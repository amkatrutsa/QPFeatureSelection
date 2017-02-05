# QPFeatureSelection

## About 

This repository presents quadratic programming feature selection approach to solve multicollinearity problem [[A. Katrutsa, V. Srijov, 2017]](http://www.sciencedirect.com/science/article/pii/S0957417417300635).
Multicollinearity in data leads to unstability, redundancy and excessive complexity of the built regression model.
To treat multicollinearity problem we state the quadratic optimization problem which is convex and has the single optimum.
This problem statement formalizes non redundancy and relevance requirements.
Thus we obtain the feature subset, which has the minimum number of similarity features and maximum number of relevant features, without any heuristics or greedy searches.

## Using in MATLAB

To test quadratic programming feature selection on test data sets run the following command in directory you choose:  
```
git clone https://github.com/amkatrutsa/QPFeatureSelection.git
cd ./mcode
```
Then in MATLAB command line from the directory `mcode` run the command:
```
main
```
After that, you will have a weight vector `x`, thresholds vector `threshold` and vectors of evaluation criteria for every value of threshold.
Now you can choose appropriate threshold and solve regression problem by standard methods which automatically give you more stable and less redundant model. 
