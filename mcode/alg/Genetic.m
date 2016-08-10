function [w] = Genetic(features, target, Params)
% Genetic algorithm.
%
% Input:
% features - matrix of features, where rows are objects, and colums are feature vectors
% target   - target feature vector
% Params   - structure with fields 
% Genetic - params for genetic; structure with fields
% nGenerations - number of generations
% nIndividuals - number of individuls of every generation
% mutationProb - probability of changing every bit of features mask
% Output
%   featuresRating - structure with rating for all features; has fields
%     isInformative - array of marks is particular feature informative (1) or not (0)
%     weight        - weight of particular feature id it is informative

nFeatures = size(features, 2); % number of features

% Default params
maxNFeatures = nFeatures; % max features number to scan
nGenerations = 10; % number of generations
nIndividuals = 10; % number of individuals of every generation
mutationProb = 0.5; % probability of changing every bit of features mask

if isfield(Params, 'Genetic') % if params for genetic exist
    if isfield(Params.Genetic, 'nGenerations') % number of generations
        nGenerations = Params.Genetic.nGenerations;
    end
    if isfield(Params.Genetic, 'nIndividuals') % number of individuls of every generation
        nIndividuals = Params.Genetic.nIndividuals;
    end
    if isfield(Params.Genetic, 'mutationProb') % probability of changing every bit of features mask
        mutationProb = Params.Genetic.mutationProb;
    end
end

generation = rand(nFeatures, nIndividuals) > 0.5; % first random generation
q = TestMask(generation, features, target, Params); % error functionals for generation

for iGeneration = 1:nGenerations
    fprintf('Generation = %d/%d\n', iGeneration, nGenerations);
    newGeneration = []; % new generation
    for iInd1 = 1:nIndividuals
        ind1 = generation(:, iInd1); % first individual
        for ind2 = generation(:, (iInd1 + 1):nIndividuals); % second individual
            partOfInd1 = rand(nFeatures, 1); % part of first individual in crossbreeding
            partOfInd2 = 1 - partOfInd1; % part of second individual in crossbreeding
            newGeneration(:, end + 1) = ind1 .* partOfInd1 + ind2 .* partOfInd2 >= 0.5; % add crossbreeding of individuals
        end
        newGeneration(:, end + 1) = xor(ind1, rand(nFeatures, 1) < mutationProb); % add mutated first individual
    end
    newGeneration = logical(newGeneration(:, FilterMasks(newGeneration, maxNFeatures))); % filter generation
    newQ = TestMask(newGeneration, features, target, Params); % error functionals for generation
    generation = [generation, newGeneration]; % concat new and olg generations
    q = [q, newQ]; % conct functionals
    idx = FilterMasks(generation, maxNFeatures); % filter it
    generation = generation(:, idx); % filtered generation
    q = q(idx); % filtered functionals 

    [~, idx] = sort(q); % sort by functional
    generation = generation(:, idx(1:nIndividuals)); % best individuals
    q = q(idx(1:nIndividuals)); % functionals for best individuals
end

featuresRating.isInformative = generation(:, 1); % informative mask
featuresRating.weights = zeros(nFeatures, 1); % start weights with zeros
featuresRating.weights(featuresRating.isInformative) = ...
    lscov(features(:, featuresRating.isInformative), target); % calculate informative weights
w = featuresRating.weights;
end