# ASOptimizer
> **ASOptimizer consists of a database and two computational models: *sequence engineering* and *chemical engineering***.

- To access ASOptimizer, we have developed a web server that runs ASOptimizer on the backend and provided free access to it. You can access it through the following link: http://asoptimizer.s-core.ai/
- To access the figures mentioned in the paper, you can use the following code
## 0. Requirements
- Python 3.6
- Install required libraries using the following command:

> <pre><code>pip3 install -r requirements.txt</code></pre>

- Please download the following data files and place npy files into the "/sequence_engineering/features" folder
  - https://drive.google.com/file/d/1RSdweiYhFJq6gs3M4ns2tH6OXzehAZZb/view?usp=sharing


## 1. Database

### 1.1 For sequence engineering
- **/ido-patent/EFO21_experiments.csv**, **/ido-patent/SK0V3_experiments.csv**
  - Database of experimental observations from the granted patent "Immunosuppression-reverting oligonucleotides inhibiting the expression of IDO"
  - We would like to express our gratitude to all contributors who provided valuable insights and resources for this research.

- **/patent_experiments/aso_features_patent.csv**
  - ASO candidates from the granted patent "Immunosuppression-reverting oligonucleotides inhibiting the expression of IDO" and corresponding features"
    
- **/features/IDO1_features.csv**
  - 19-mer ASO candidates for regulating IDO1 mRNA and corresponding features
   

## 2. Sequence Engineering

To perform sequence engineering, use the following command:
> <pre><code>python3 main.py --mode 'train' --target 'ido1' --seq_len 19 --num_candidates 6 --rnastructure 'mfold' </code></pre>

To perform sequence engineering with *your chosen setting*, use the following command:
> <pre><code>python3 main.py --mode 'eval' --target 'ido1' --seq_len 19 --num_candidates 6 --rnastructure 'mfold' --a_star "1 1 1 1" </code></pre>

### Code

- **Linear_regression.ipynb**: 
  - Displays a scatter plot comparing experimentally observed inhibition rates (x-axis) with their predicted values (y-axis) (Figure 1A and Figure 1B)
  - Comparison with sfold

- **Plot_Contour.ipynb**: 
  - Displays a surface plot of Pearson correlation (ρ), represented by contour plots (Figure 1C and Figure 1D)

- **Getting_ASOpt_top6.ipynb**: 
  - Shows a histogram depicting the predicted scores of complementary ASOs, each 19 nucleotides in length, targeting the IDO1 gene (Figure S2)
