% Input IDX_MATRIX: 
%   rows number is repeat times,
%   columns number is items number of basic data
%
function [MRNA, PROTEIN, PROTEIN_ORIGINAL, T, N, names] =  GeneDataLoad()
load hmec;
MRNA = data{1};
P1 = data{2};
MRNA = normalize(MRNA)';
PROTEIN = normalize(P1(:,2:7))';
PROTEIN_ORIGINAL = normalize(P1)';
[T,N] = size(MRNA);
names = names;