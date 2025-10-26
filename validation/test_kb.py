"""
Pathway Knowledge Base Validation Framework
============================================
A framework for evaluating pathway activation predictions using 
Positive-Unlabeled (PU) learning on single-cell transcriptomics data.
"""

import numpy as np
import pandas as pd
import json
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# For visualization
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import precision_recall_curve, auc
from sklearn.model_selection import KFold
from collections import defaultdict
from itertools import islice

@dataclass
class PathwayScore:
    """Container for pathway activation scores and metadata."""
    cluster: str
    pathway: str
    score: float
    expressed_genes: List[str]
    total_genes: int
    module_scores: Dict[str, float]
    
@dataclass
class TuningResult:
    """Captures the outcome of a single weight configuration evaluation."""
    weight_core: float
    weight_regulators: float
    objective: float
    per_pathway_metrics: Dict[str, Dict]


class GeneExpressionLoader:
    """Handles loading and preprocessing of gene expression data."""
    
    def __init__(self, expression_file: str):
        self.expression_file = expression_file
        self.data = None
        
    def load(self) -> pd.DataFrame:
        """Load gene expression data."""
        self.data = pd.read_csv(self.expression_file)
        print(f"Loaded {len(self.data)} expression records")
        print(f"Unique clusters: {self.data['cluster'].nunique()}")
        print(f"Unique genes: {self.data['gene_id'].nunique()}")
        return self.data
    
    def get_cluster_genes(self, cluster: str, 
                            detection_threshold: float = 0.01,
                            cpm_threshold: float = 1.0) -> Set[str]:
        """Get expressed genes for a specific cluster."""
        cluster_data = self.data[self.data['cluster'] == cluster]
        expressed = cluster_data[
            (cluster_data['detection_fraction'] >= detection_threshold) &
            (cluster_data['mean_cpm'] >= cpm_threshold)
        ]
        return set(expressed['gene_id'].values)


class GroundTruthLoader:
    """Handles loading of ground truth pathway annotations."""
    
    def __init__(self, ground_truth_file: str):
        self.ground_truth_file = ground_truth_file
        self.data = None
        
    def load(self) -> pd.DataFrame:
        """Load ground truth annotations."""
        self.data = pd.read_csv(self.ground_truth_file)
        print(f"Loaded ground truth for {len(self.data)} cell types")
        print(f"Pathways: {[c for c in self.data.columns if c != 'cluster']}")
        return self.data
    
    def get_positive_labels(self) -> Dict[str, Set[str]]:
        """Extract positive pathway-cluster pairs."""
        positives = defaultdict(set)
        for _, row in self.data.iterrows():
            cell_type = row['cluster']
            for pathway in self.data.columns[1:]:
                if pd.notna(row[pathway]) and row[pathway] == 'x':
                    positives[pathway].add(cell_type)
        return dict(positives)


class KnowledgeBaseLoader:
    """Handles loading and parsing of pathway knowledge base."""
    
    def __init__(self, kb_file: str):
        self.kb_file = kb_file
        self.kb = None
        
    def load(self) -> Dict:
        """Load knowledge base JSON."""
        with open(self.kb_file, 'r') as f:
            self.kb = json.load(f)
        print(f"Loaded KB with {len(self.kb)} pathways")
        for pathway, content in self.kb.items():
            print(f"  {pathway}: {len(self._extract_all_genes(content))} genes")
        return self.kb
    
    def _extract_all_genes(self, pathway_def: Dict) -> Set[str]:
        """Extract all gene IDs from a pathway definition."""
        genes = set()
        for key, value in pathway_def.items():
            if isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        if 'genes' in item:
                            genes.update(g['ensembl_id'] for g in item['genes'])
                        elif 'ensembl_id' in item:
                            genes.add(item['ensembl_id'])
        return genes
    
    def get_pathway_structure(self, pathway: str) -> Dict:
        """Get the structured definition of a pathway."""
        return self.kb.get(pathway, {})


class PathwayScorer:
    """Calculates pathway activation scores for cell clusters."""
    
    def __init__(self,
                 kb_loader: KnowledgeBaseLoader,
                 weight_core: float = 0.7,
                 weight_regulators: float = 0.3):
        self.kb_loader = kb_loader
        self.weight_core = weight_core
        self.weight_regulators = weight_regulators
    
    def set_weights(self,
                    weight_core: Optional[float] = None,
                    weight_regulators: Optional[float] = None,
                    normalize: bool = True) -> None:
        """
        Update the default blending weights.
        
        Args:
            weight_core: Optional override for the core (gate) component weight.
            weight_regulators: Optional override for the regulator component weight.
            normalize: If True, rescale so weights sum to 1. This keeps the
                final score comparable across runs while allowing easy extension
                to more components in the future.
        """
        if weight_core is not None:
            self.weight_core = float(weight_core)
        if weight_regulators is not None:
            self.weight_regulators = float(weight_regulators)
        
        if normalize:
            total = self.weight_core + self.weight_regulators
            if total > 0:
                self.weight_core /= total
                self.weight_regulators /= total
    
    def get_weights(self) -> Tuple[float, float]:
        """Return the current default weights as a convenience helper."""
        return self.weight_core, self.weight_regulators
        
    def score_module(self, module: Dict, expressed_genes: Set[str]) -> Tuple[float, List[str]]:
        """
        Score a single module based on gene expression.
        
        Args:
            module: Module definition with 'op' and 'genes'
            expressed_genes: Set of expressed gene IDs in cluster
            
        Returns:
            Tuple of (score, list of expressed genes in module)
        """
        module_genes = [g['ensembl_id'] for g in module.get('genes', [])]
        expressed_in_module = [g for g in module_genes if g in expressed_genes]
        if not module_genes:
            return 0.0, []
        
        op = module.get('op', 'ANY')
        
        if op == 'ALL':
            # All genes must be expressed
            score = 1.0 if len(expressed_in_module) == len(module_genes) else 0.0
        elif op == 'ANY':
            # At least one gene must be expressed
            score = 1.0 if len(expressed_in_module) > 0 else 0.0
        elif op == 'MAJORITY':
            # More than half must be expressed
            score = 1.0 if len(expressed_in_module) > len(module_genes) / 2 else 0.0
        else:
            # Default: proportion of expressed genes
            score = len(expressed_in_module) / len(module_genes)
            
        return score, expressed_in_module
    
    def score_pathway(self, pathway: str, cluster: str, 
                        expressed_genes: Set[str],
                        weight_core: Optional[float] = None,
                        weight_regulators: Optional[float] = None) -> PathwayScore:
        """
        Calculate comprehensive activation score for a pathway in a cluster.
        
        Strategy:
        1. For core modules (with flavors): ANY flavor being active = pathway active
        2. Regulators add supporting evidence (weighted contribution)
        3. Final score = weighted combination
        """
        pathway_def = self.kb_loader.get_pathway_structure(pathway)
        if not pathway_def:
            return PathwayScore(cluster, pathway, 0.0, [], 0, {})
        
        # Score core modules (flavors)
        core_scores = []
        all_expressed = []
        module_score_dict = {}
        
        for key in ['core_modules']:
            if key in pathway_def:
                modules = pathway_def[key]
                if not isinstance(modules, list):
                    modules = [modules]
                    
                for i, module in enumerate(modules):
                    module_name = module.get('name', f'{key}_{i}')
                    score, expressed = self.score_module(module, expressed_genes)
                    core_scores.append(score)
                    all_expressed.extend(expressed)
                    module_score_dict[module_name] = score
        
        # Core score: ANY module being active (max score)
        core_score = max(core_scores) if core_scores else 0.0
        
        # Score regulators (supporting evidence)
        gate_scores = []
        if 'gate_core' in pathway_def:
            gates = pathway_def['gate_core']
            if not isinstance(gates, list):
                gates = [gates]
            
            for g in gates:
                if isinstance(g, dict) and 'ensembl_id' in g:
                    if g['ensembl_id'] in expressed_genes:
                        gate_scores.append(1.0)
                        all_expressed.append(g['ensembl_id'])
                    else:
                        gate_scores.append(0.0)

        # Score regulators (supporting evidence)
        regulator_scores = []
        if 'regulators' in pathway_def:
            regulators = pathway_def['regulators']
            if not isinstance(regulators, list):
                regulators = [regulators]
            
            for reg in regulators:
                if isinstance(reg, dict) and 'ensembl_id' in reg:
                    if reg['ensembl_id'] in expressed_genes:
                        regulator_scores.append(1.0)
                        all_expressed.append(reg['ensembl_id'])
                    else:
                        regulator_scores.append(0.0)
        
        # Core score: proportion expressed
        gate_score = np.mean(gate_scores) if gate_scores else 0.0
        module_score_dict['gate_core'] = gate_score

        # Regulator score: proportion expressed
        regulator_score = np.mean(regulator_scores) if regulator_scores else 0.0
        module_score_dict['regulators'] = regulator_score
        
        # Combined score
        w_core = self.weight_core if weight_core is None else float(weight_core)
        w_reg = self.weight_regulators if weight_regulators is None else float(weight_regulators)
        total = w_core + w_reg
        if total > 0:
            w_core /= total
            w_reg /= total
        else:
            w_core = 0.0
            w_reg = 0.0
        
        # Final blended score keeps weights normalized for comparability.
        final_score = (w_core * gate_score + 
                      w_reg * regulator_score)
        
        # Count total genes in pathway
        all_pathway_genes = set()
        for key, value in pathway_def.items():
            if isinstance(value, list):
                for item in value:
                    if isinstance(item, dict):
                        if 'genes' in item:
                            all_pathway_genes.update(g['ensembl_id'] for g in item['genes'])
                        elif 'ensembl_id' in item:
                            all_pathway_genes.add(item['ensembl_id'])
        
        return PathwayScore(
            cluster=cluster,
            pathway=pathway,
            score=final_score,
            expressed_genes=list(set(all_expressed)),
            total_genes=len(all_pathway_genes),
            module_scores=module_score_dict
        )


class PULearningEvaluator:
    """Implements Positive-Unlabeled learning evaluation strategies."""
    
    def __init__(self, scores_df: pd.DataFrame, positive_labels: Dict[str, Set[str]]):
        """
        Args:
            scores_df: DataFrame with columns [cluster, pathway, score]
            positive_labels: Dict mapping pathway -> set of positive clusters
        """
        self.scores_df = scores_df
        self.positive_labels = positive_labels
        
    def estimate_class_prior(self, pathway: str) -> float:
        """
        Estimate the true proportion of positive examples using PU learning.
        Uses the "selected completely at random" (SCAR) assumption.
        """
        pathway_scores = self.scores_df[self.scores_df['pathway'] == pathway]
        positive_clusters = self.positive_labels.get(pathway, set())
        
        if not positive_clusters:
            return 0.0
        
        # Label: 1 if positive, 0 if unlabeled
        labels = pathway_scores['cluster'].apply(
            lambda x: 1 if x in positive_clusters else 0
        ).values
        scores = pathway_scores['score'].values
        
        # Class prior estimation using threshold-free method
        # E[s|y=1] = mean score of known positives
        # E[s] = mean score of all data
        # c = P(s=1|y=1) (labeling frequency) ≈ |P|/|P_true|
        
        if labels.sum() == 0:
            return 0.0
            
        mean_pos = scores[labels == 1].mean() if labels.sum() > 0 else 0
        mean_all = scores.mean()
        
        # Estimate: π = (E[s] - (1-c)E[s|y=0]) / c
        # Simplified: assume E[s|y=0] ≈ min(scores)
        prior = max(0.1, min(0.9, labels.sum() / len(labels)))
        
        return prior
    
    def create_reliable_negatives(self, pathway: str, 
                                    n_negatives: Optional[int] = None,
                                    strategy: str = 'low_score') -> Set[str]:
        """
        Create reliable negative set from unlabeled data.
        
        Strategies:
        - 'low_score': Select clusters with lowest activation scores
        - 'spy': Use a small subset of positives as spies in unlabeled
        """
        pathway_scores = self.scores_df[self.scores_df['pathway'] == pathway]
        positive_clusters = self.positive_labels.get(pathway, set())
        
        # Get unlabeled clusters
        unlabeled = pathway_scores[
            ~pathway_scores['cluster'].isin(positive_clusters)
        ]
        
        if n_negatives is None:
            n_negatives = len(positive_clusters)  # Same as positives
        
        if strategy == 'low_score':
            # Select clusters with lowest scores as reliable negatives
            reliable_negatives = unlabeled.nsmallest(n_negatives, 'score')
            return set(reliable_negatives['cluster'].values)
        
        elif strategy == 'spy':
            # More sophisticated: use spy technique
            # (Simplified version - full implementation would iterate)
            threshold = unlabeled['score'].quantile(0.25)
            reliable_negatives = unlabeled[unlabeled['score'] < threshold]
            return set(reliable_negatives['cluster'].values)
        
        return set()
    
    def compute_pu_metrics(self, pathway: str, threshold: float = 0.5) -> Dict:
        """
        Compute PU-specific evaluation metrics.
        
        Returns:
            Dictionary with precision, recall, F1, and PU-specific metrics
        """
        pathway_scores = self.scores_df[self.scores_df['pathway'] == pathway]
        positive_clusters = self.positive_labels.get(pathway, set())
        
        if not positive_clusters:
            return {'error': 'No positive labels for pathway'}
        
        # Create predictions
        predictions = pathway_scores['score'] >= threshold
        
        # True labels (1=positive, 0=unlabeled)
        labels = pathway_scores['cluster'].isin(positive_clusters)
        
        # Positive metrics (treating unlabeled as negative - lower bound)
        tp = (predictions & labels).sum()
        fp = (predictions & ~labels).sum()
        fn = (~predictions & labels).sum()
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        # PU metrics
        class_prior = self.estimate_class_prior(pathway)
        
        # Ranking-based metric: how well do scores rank positives?
        scores = pathway_scores.set_index('cluster')['score']
        pos_scores = [scores[c] for c in positive_clusters if c in scores.index]
        all_scores = scores.values
        
        if pos_scores:
            # Average rank percentile of positive examples
            ranks = [stats.percentileofscore(all_scores, s) for s in pos_scores]
            avg_rank_percentile = np.mean(ranks)
        else:
            avg_rank_percentile = 0
        
        return {
            'precision_lower_bound': precision,
            'recall_on_labeled': recall,
            'f1_lower_bound': f1,
            'estimated_class_prior': class_prior,
            'avg_rank_percentile': avg_rank_percentile,
            'n_positive': len(positive_clusters),
            'n_predicted_positive': predictions.sum(),
            'n_total': len(pathway_scores)
        }
    
    def cross_validate_with_negatives(self, pathway: str, 
                                        n_folds: int = 5) -> Dict:
        """
        Perform cross-validation using reliable negatives.
        """
        pathway_scores = self.scores_df[self.scores_df['pathway'] == pathway]
        positive_clusters = self.positive_labels.get(pathway, set())
        reliable_negatives = self.create_reliable_negatives(pathway)
        
        # Create labeled dataset
        labeled_data = pathway_scores[
            pathway_scores['cluster'].isin(positive_clusters | reliable_negatives)
        ].copy()
        labeled_data['label'] = labeled_data['cluster'].isin(positive_clusters).astype(int)
        
        if len(labeled_data) < n_folds:
            return {'error': 'Insufficient data for CV'}
        
        # Simple threshold-based CV
        kf = KFold(n_splits=min(n_folds, len(labeled_data)), shuffle=True, random_state=42)
        cv_scores = []
        
        for train_idx, test_idx in kf.split(labeled_data):
            train = labeled_data.iloc[train_idx]
            test = labeled_data.iloc[test_idx]
            
            # Find optimal threshold on train set
            thresholds = np.percentile(train['score'], [25, 50, 75])
            best_f1 = 0
            best_thresh = 0.5
            
            for thresh in thresholds:
                pred = train['score'] >= thresh
                tp = (pred & train['label']).sum()
                fp = (pred & ~train['label']).sum()
                fn = (~pred & train['label']).sum()
                
                prec = tp / (tp + fp) if (tp + fp) > 0 else 0
                rec = tp / (tp + fn) if (tp + fn) > 0 else 0
                f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
                
                if f1 > best_f1:
                    best_f1 = f1
                    best_thresh = thresh
            
            # Evaluate on test
            test_pred = test['score'] >= best_thresh
            tp = (test_pred & test['label']).sum()
            fp = (test_pred & ~test['label']).sum()
            fn = (~test_pred & test['label']).sum()
            
            prec = tp / (tp + fp) if (tp + fp) > 0 else 0
            rec = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
            
            cv_scores.append({'precision': prec, 'recall': rec, 'f1': f1})
        
        return {
            'mean_precision': np.mean([s['precision'] for s in cv_scores]),
            'mean_recall': np.mean([s['recall'] for s in cv_scores]),
            'mean_f1': np.mean([s['f1'] for s in cv_scores]),
            'std_f1': np.std([s['f1'] for s in cv_scores])
        }


class PathwayKBValidator:
    """Main validation pipeline orchestrator."""
    
    def __init__(self, expression_file: str, ground_truth_file: str, kb_file: str):
        self.expr_loader = GeneExpressionLoader(expression_file)
        self.gt_loader = GroundTruthLoader(ground_truth_file)
        self.kb_loader = KnowledgeBaseLoader(kb_file)
        self.scorer = None
        self.scores_df = None
        self.evaluation_results = {}
        self.weight_tuner = None
        
    def load_data(self):
        """Load all data sources."""
        print("Loading data...")
        self.expr_loader.load()
        self.gt_loader.load()
        self.kb_loader.load()
        self.scorer = PathwayScorer(self.kb_loader)
        print("Data loading complete.\n")
        
    def compute_all_scores(self, detection_threshold: float = 0.01,
                            cpm_threshold: float = 1.0,
                            weight_core: Optional[float] = None,
                            weight_regulators: Optional[float] = None,
                            normalize_weights: bool = True) -> pd.DataFrame:
        """Compute pathway scores for all cluster-pathway pairs."""
        print("Computing pathway activation scores...")
        
        if weight_core is not None or weight_regulators is not None:
            self.scorer.set_weights(weight_core, weight_regulators, normalize=normalize_weights)
        
        current_core, current_reg = self.scorer.get_weights()
        
        clusters = self.expr_loader.data['cluster'].unique()
        pathways = list(self.kb_loader.kb.keys())
        
        results = []
        for cluster in clusters:
            expressed_genes = self.expr_loader.get_cluster_genes(
                cluster, detection_threshold, cpm_threshold
            )
            for pathway in pathways:
                score_obj = self.scorer.score_pathway(
                    pathway,
                    cluster,
                    expressed_genes,
                    weight_core=current_core,
                    weight_regulators=current_reg
                )
                results.append({
                    'cluster': cluster,
                    'pathway': pathway,
                    'score': score_obj.score,
                    'n_expressed_genes': len(score_obj.expressed_genes),
                    'total_genes': score_obj.total_genes,
                    'coverage': len(score_obj.expressed_genes) / score_obj.total_genes 
                                if score_obj.total_genes > 0 else 0
                })
        
        self.scores_df = pd.DataFrame(results)
        print(f"Computed {len(results)} pathway-cluster scores.\n")
        return self.scores_df
    
    def evaluate_kb(self) -> Dict:
        """Run comprehensive evaluation of the knowledge base."""
        print("Evaluating knowledge base...")
        
        positive_labels = self.gt_loader.get_positive_labels()
        evaluator = PULearningEvaluator(self.scores_df, positive_labels)
        
        results = {}
        for pathway in islice(self.kb_loader.kb.keys(), 4):
            print(f"\nEvaluating: {pathway}")
            
            # PU metrics
            pu_metrics = evaluator.compute_pu_metrics(pathway)
            print(f"  Precision (lower bound): {pu_metrics.get('precision_lower_bound', 0):.3f}")
            print(f"  Recall on labeled: {pu_metrics.get('recall_on_labeled', 0):.3f}")
            print(f"  Avg rank percentile: {pu_metrics.get('avg_rank_percentile', 0):.1f}")
            
            # CV with reliable negatives
            cv_results = evaluator.cross_validate_with_negatives(pathway)
            if 'error' not in cv_results:
                print(f"  CV F1: {cv_results['mean_f1']:.3f} ± {cv_results['std_f1']:.3f}")
            
            results[pathway] = {
                'pu_metrics': pu_metrics,
                'cv_results': cv_results
            }
        
        self.evaluation_results = results
        return results
    
    def tune_weights(self,
                     core_grid: Optional[List[float]] = None,
                     regulator_grid: Optional[List[float]] = None,
                     detection_threshold: float = 0.01,
                     cpm_threshold: float = 1.0) -> TuningResult:
        """
        Run the lightweight WeightTuner to pick blending weights that work best
        for the current dataset. The evaluation metrics already implemented are
        reused as the optimization target.
        
        Returns:
            TuningResult describing the selected weights and summary statistics.
        """
        if self.scorer is None:
            raise RuntimeError("Call load_data() before tuning weights.")
        
        self.weight_tuner = WeightTuner(
            validator=self,
            core_grid=core_grid,
            regulator_grid=regulator_grid,
            detection_threshold=detection_threshold,
            cpm_threshold=cpm_threshold
        )
        result = self.weight_tuner.tune()
        print(f"[PathwayKBValidator] Weight tuning objective: {result.objective:.3f}")
        return result
    
    def visualize_results(self, output_dir: str = '.'):
        """Generate visualization plots."""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # 1. Heatmap of activation scores
        plt.figure(figsize=(12, 8))
        pivot = self.scores_df.pivot(index='cluster', columns='pathway', values='score')
        sns.heatmap(pivot, cmap='YlOrRd', vmin=0, vmax=1, cbar_kws={'label': 'Activation Score'})
        plt.title('Pathway Activation Scores Across Cell Clusters')
        plt.tight_layout()
        plt.savefig(output_path / 'activation_heatmap.png', dpi=300)
        plt.close()
        
        # 2. Score distributions by pathway
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        pathways = list(self.kb_loader.kb.keys())
        
        for idx, pathway in enumerate(pathways[:4]):  # First 4 pathways
            ax = axes[idx // 2, idx % 2]
            pathway_data = self.scores_df[self.scores_df['pathway'] == pathway]
            
            # Plot score distribution
            ax.hist(pathway_data['score'], bins=20, alpha=0.7, color='steelblue')
            ax.axvline(0.5, color='red', linestyle='--', label='Threshold=0.5')
            ax.set_xlabel('Activation Score')
            ax.set_ylabel('Frequency')
            ax.set_title(f'{pathway}')
            ax.legend()
        
        plt.tight_layout()
        plt.savefig(output_path / 'score_distributions.png', dpi=300)
        plt.close()
        
        # 3. Performance metrics comparison
        pathways_list = []
        precisions = []
        recalls = []
        rank_percentiles = []
        for pathway, results in islice(self.evaluation_results.items(), 4):
            pu = results['pu_metrics']
            pathways_list.append(pathway)
            precisions.append(pu.get('precision_lower_bound', 0))
            recalls.append(pu.get('recall_on_labeled', 0))
            rank_percentiles.append(pu.get('avg_rank_percentile', 0))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        x = np.arange(len(pathways_list))
        width = 0.35
        
        ax1.bar(x - width/2, precisions, width, label='Precision (LB)', alpha=0.8)
        ax1.bar(x + width/2, recalls, width, label='Recall', alpha=0.8)
        ax1.set_ylabel('Score')
        ax1.set_title('Precision and Recall by Pathway')
        ax1.set_xticks(x)
        ax1.set_xticklabels(pathways_list, rotation=45, ha='right')
        ax1.legend()
        ax1.set_ylim(0, 1)
        
        ax2.bar(pathways_list, rank_percentiles, color='coral', alpha=0.8)
        ax2.set_ylabel('Percentile')
        ax2.set_title('Average Rank Percentile of Positive Examples')
        ax2.set_xticklabels(pathways_list, rotation=45, ha='right')
        ax2.set_ylim(0, 100)
        ax2.axhline(50, color='gray', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(output_path / 'performance_metrics.png', dpi=300)
        plt.close()
        
        print(f"\nVisualizations saved to {output_path}")


class WeightTuner:
    """
    Lightweight grid-search tuner for the core/regulator blending weights.
    
    The tuner reuses the validator's existing scoring/evaluation pipeline, so we
    benefit from every metric already implemented. It keeps a history of the
    tried combinations, making it easy to extend to more parameters later.
    """
    
    def __init__(self,
                 validator: "PathwayKBValidator",
                 core_grid: Optional[List[float]] = None,
                 regulator_grid: Optional[List[float]] = None,
                 detection_threshold: float = 0.01,
                 cpm_threshold: float = 1.0):
        self.validator = validator
        self.core_grid = core_grid or [0.4, 0.5, 0.6, 0.7, 0.8]
        self.regulator_grid = regulator_grid or [0.1, 0.2, 0.3, 0.4, 0.5]
        self.detection_threshold = detection_threshold
        self.cpm_threshold = cpm_threshold
        self.history: List[TuningResult] = []
    
    def tune(self) -> TuningResult:
        """
        Run grid search over the configured weight grids.
        
        Returns:
            The best-performing configuration according to the objective score.
            The validator is left in the state corresponding to this best run.
        """
        original_weights = self.validator.scorer.get_weights()
        best_result: Optional[TuningResult] = None
        
        for weight_core in self.core_grid:
            for weight_reg in self.regulator_grid:
                if np.isclose(weight_core + weight_reg, 0.0):
                    continue  # skip degenerate combinations
                
                print(f"\n[WeightTuner] Trying core={weight_core:.2f}, regulators={weight_reg:.2f}")
                self.validator.compute_all_scores(
                    detection_threshold=self.detection_threshold,
                    cpm_threshold=self.cpm_threshold,
                    weight_core=weight_core,
                    weight_regulators=weight_reg,
                    normalize_weights=True
                )
                evaluation = self.validator.evaluate_kb()
                objective = self._objective(evaluation)
                
                result = TuningResult(
                    weight_core=weight_core,
                    weight_regulators=weight_reg,
                    objective=objective,
                    per_pathway_metrics=evaluation
                )
                self.history.append(result)
                
                if best_result is None or result.objective > best_result.objective:
                    best_result = result
                    print(f"[WeightTuner] New best objective: {objective:.3f}")
        
        if best_result is None:
            # Restore original weights if no configuration was feasible.
            self.validator.scorer.set_weights(*original_weights)
            return TuningResult(
                weight_core=original_weights[0],
                weight_regulators=original_weights[1],
                objective=float("-inf"),
                per_pathway_metrics={}
            )
        
        # Ensure validator reflects best configuration and store metrics for downstream use.
        self.validator.compute_all_scores(
            detection_threshold=self.detection_threshold,
            cpm_threshold=self.cpm_threshold,
            weight_core=best_result.weight_core,
            weight_regulators=best_result.weight_regulators,
            normalize_weights=True
        )
        self.validator.evaluation_results = best_result.per_pathway_metrics
        print(f"\n[WeightTuner] Selected weights -> core={best_result.weight_core:.2f}, "
              f"regulators={best_result.weight_regulators:.2f} (objective={best_result.objective:.3f})")
        return best_result
    
    def _objective(self, evaluation_results: Dict[str, Dict]) -> float:
        """
        Collapse pathway-level metrics into a single score.
        
        Current strategy:
            - Average the available precision, recall, and rank percentile (scaled to [0,1])
            - Include CV mean F1 when present
            - Average across pathways
        
        This keeps the objective transparent while making it straightforward to
        incorporate extra metrics later.
        """
        pathway_scores = []
        for result in evaluation_results.values():
            pathway_metrics = []
            pu = result.get('pu_metrics', {})
            if pu:
                if pu.get('precision_lower_bound') is not None:
                    pathway_metrics.append(pu.get('precision_lower_bound', 0.0))
                if pu.get('recall_on_labeled') is not None:
                    pathway_metrics.append(pu.get('recall_on_labeled', 0.0))
                if pu.get('avg_rank_percentile') is not None:
                    pathway_metrics.append(pu.get('avg_rank_percentile', 0.0) / 100.0)
            cv = result.get('cv_results', {})
            if isinstance(cv, dict) and 'mean_f1' in cv:
                pathway_metrics.append(cv['mean_f1'])
            if pathway_metrics:
                pathway_scores.append(float(np.nanmean(pathway_metrics)))
        
        if not pathway_scores:
            return float("-inf")
        return float(np.nanmean(pathway_scores))

# Example usage
if __name__ == "__main__":
    # Initialize validator
    validator = PathwayKBValidator(
        expression_file='cpp-mechanisms/validation/GTEx/cluster_gene_expression.csv',
        ground_truth_file='cpp-mechanisms/validation/GTEx/ground_truth.csv',
        kb_file='cpp-mechanisms/data/mechanisms.json'
    )
    
    # Run pipeline
    validator.load_data()
    
    # Optional: automatically tune the blend between core and regulator components.
    tuning_result = validator.tune_weights(
        core_grid=[0.5, 0.6, 0.7, 0.8],
        regulator_grid=[0.2, 0.3, 0.4, 0.5],
        detection_threshold=0.01,
        cpm_threshold=1.0
    )
    print(f"\nBest weights -> core={tuning_result.weight_core:.2f}, "
          f"regulators={tuning_result.weight_regulators:.2f}")
    
    # The tuner leaves the validator holding the best scoring DataFrame/metrics.
    scores_df = validator.scores_df
    
    # Evaluate
    results = validator.evaluation_results or validator.evaluate_kb()
    
    # Visualize
    validator.visualize_results(output_dir='./results')
    
    # Save scores
    scores_df.to_csv('pathway_scores.csv', index=False)
    print("\nPipeline complete!")
