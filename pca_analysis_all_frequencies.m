% PCA analysis of all fly results irrespective of frequency

close all; clear;

load results

load slrp_lrpr.mat

struct_names = {'freq1dot25hz','freq6dot25hz','freq12dot5hz','freq25hz'};

X = [];

%concatenate data
for f = 1:4
   
    data = results.(struct_names{f});
    
    for i = 1:length(data)
       
        X = [X;normalize(data(i).seq_eff_profile).']; %#ok<AGROW>
%           X = [X;normalize(rand(1,16))]; %#ok<AGROW>
    end

end

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(X);

target = [normalize(slrp).'; normalize(-lrpr).'; normalize(weird).'];
% target = [rand(1,16); rand(1,16); rand(1,16)];

B = rotatefactors(COEFF(1:3,:),'Method','procrustes','Target',target,'type','oblique');

for i = 1:3
   figure; create_seq_eff_plot(real(B(i,:)).',[]);
end