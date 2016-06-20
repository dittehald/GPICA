function defaultsetonoffCONV

global h_main

if get(h_main.settings.ica.checkconv,'Value') == 1  % Perform PCA on
    set(h_main.settings.ica.algorithm.convlags,'Enable','on')
%     set(h_main.settings.ica.algorithm.convlagstext,'Enable','on')
    
elseif get(h_main.settings.ica.checkconv,'Value') == 0
       set(h_main.settings.ica.algorithm.convlags,'Enable','off')
%        set(h_main.settings.ica.algorithm.convlagstext,'Enable','off')
end