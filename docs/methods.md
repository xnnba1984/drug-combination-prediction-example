# Analysis Methods

## Prediction Task

One observation is an unordered drug pair in a cell line. The source labels contain 23,052 rows. Averaging the 315 repeated triplets gives 22,737 observations: 583 pairs across 39 cell lines. Each drug appears in other training pairs, and every training fold includes all 39 cell lines. The prediction task concerns new pairings of these known drugs in this panel.

The endpoint is the continuous Loewe synergy score in the public DeepSynergy labels. CSS is joined for retrospective interpretation after ranking. It is excluded from feature construction, model fitting, tuning, and candidate selection.

## Predictors

Drugs have 4,387 published features, and cell lines have 3,984 gene expression features. Constant columns are removed. Principal component analysis (PCA) is fitted separately to the fixed inventory of 38 drug profiles and 39 cell line profiles, retaining at least 80% of variance. This gives 23 drug and 14 cell line components. All these entity profiles are available in the stated deployment task. PCA uses no outcome measurements.

For a pair, the drug components enter as means and absolute differences. This makes predictions invariant to exchanging Drug 1 and Drug 2. Adding the cell line components gives 60 base predictors. The linear feature models also include products of the 46 pair components and 14 cell line components, for 704 predictors in total.

The identity design has 1,559 columns: 38 drug indicators, 39 cell line indicators, and 1,482 drug-cell line interactions. Each observation has indicators for both drugs. The hybrid design combines the identity terms with the 704 feature terms, giving 2,263 columns.

## Validation And Models

The analysis retains the five published pair-level folds. One fold is held out for evaluation at a time. Tuning uses only the other four folds, keeping each pair's 39 observations together in the inner splits as well.

| Model | Specification |
| :--- | :--- |
| Training mean | Mean outcome in the outer training set |
| Identity ridge | Identity design; `glmnet`, alpha 0 |
| Feature elastic net | Linear feature design; alpha selected from 0, 0.5, and 1 |
| Hybrid elastic net | Combined identity and feature design; the same alpha grid |
| Feature random forest | Base features; `ranger`, mtry 8 or 20 and minimum node size 5 or 20 |

The linear models use 80 candidate lambda values and choose `lambda.min` by inner-fold mean squared error. `glmnet` standardizes its predictors during fitting. Random forest tuning uses 250 trees per fit; the chosen setting is refitted with 750 trees. Explicit forest seeds are determined by the outer fold, tuning-grid entry, and inner fold. Bootstrap and permutation checks use seed 20260809.

## Evaluation And Ranking

RMSE is calculated across all 22,737 out-of-fold predictions. Supporting measures are mean absolute error and Pearson and Spearman correlations. The 95% intervals use 1,000 bootstrap samples of drug pairs, with each sampled pair contributing all 39 observations. These intervals summarize variation in the fixed out-of-fold predictions; model fitting is not repeated inside the bootstrap.

Triplets are ranked by prediction within each outer fold. The top 10% is rounded up to a whole number of observations. A Loewe score above 10 is used to summarize the observed positive rate. Enrichment is calculated separately in each fold and averaged with weights proportional to the selected count. It therefore differs from dividing the two pooled percentages.

For pair-level ranking, the three cell lines with the highest predicted scores are selected for every held-out pair. Their mean prediction is the pair score. The top 10% of pairs in each fold, again rounded up, gives 60 displayed pairs. Observed Loewe scores and CSS values are averaged over those same three cell lines after selection. The reported mean difference of 8.55 uses all 60 pairs, including those with observed values below prediction.

Twenty negative-control runs permute complete 39-cell-line outcome profiles across pairs. Each run repeats identity-ridge fitting and tuning with 60 candidate lambda values. This preserves the cell line structure while disrupting the link between pair identity and the outcome profile.

## Interpretation

The identity ridge model performs best for this screen and validation target. Its candidate list is concentrated around a few drugs. The comparison supports using an identity baseline when evaluating molecular features for new combinations of known drugs. Independent experiments are needed to assess reproducibility and biological relevance of the prioritized pairs.

The analysis does not evaluate generalization to an entirely unseen drug, a new cell line, or another screen. Those questions require corresponding validation splits and training-only feature preparation.
