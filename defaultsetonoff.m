function defaultsetonoff

global h_main ica

if get(h_main.settings.ica.check,'Value') == 0  % Perform PCA on
    set(h_main.settings.ica.algorithm.ss,'Enable','on')
    set(h_main.settings.ica.algorithm.sr,'Enable','on')
    set(h_main.settings.ica.algorithm.ts,'Enable','on')
    set(h_main.settings.ica.algorithm.tr,'Enable','on')
    set(h_main.settings.ica.algorithm.emean,'Enable','on')
%     set(h_main.settings.ica.algorithm.er,'Enable','on')
    set(h_main.settings.ica.algorithm.ess,'Enable','on')
     set(h_main.settings.ica.algorithm.std,'Enable','on')
elseif get(h_main.settings.ica.check,'Value') == 1
       set(h_main.settings.ica.algorithm.ss,'Enable','off')
    set(h_main.settings.ica.algorithm.sr,'Enable','off')
    set(h_main.settings.ica.algorithm.ts,'Enable','off')
    set(h_main.settings.ica.algorithm.tr,'Enable','off')
    set(h_main.settings.ica.algorithm.emean,'Enable','off')
   % set(h_main.settings.ica.algorithm.er,'Enable','off')
    set(h_main.settings.ica.algorithm.ess,'Enable','off')
     set(h_main.settings.ica.algorithm.std,'Enable','on')
    
end