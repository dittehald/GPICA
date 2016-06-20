global h_main ica
temp = str2num(ica.TR{1})*1.2;

set(h_main.settings.ica.algorithm.emean,'string',temp);