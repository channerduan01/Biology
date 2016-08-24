%%
% Data Generator
%

function [MRNA, PROTEIN] = AddNoise(MRNA_ORIGINAL, PROTEIN_ORIGINAL, amplitude)
MRNA = MRNA_ORIGINAL;
PROTEIN = PROTEIN_ORIGINAL;

% Add noise
NOISE_1 = randn(size(MRNA))*amplitude;
NOISE_2 = randn(size(PROTEIN))*amplitude;
MRNA(MRNA>0) = MRNA(MRNA>0) + NOISE_1(MRNA>0);
PROTEIN(PROTEIN>0) = PROTEIN(PROTEIN>0) + NOISE_2(PROTEIN>0);

% Normalize original data
MRNA = normalize(MRNA);
PROTEIN = normalize(PROTEIN);
end




