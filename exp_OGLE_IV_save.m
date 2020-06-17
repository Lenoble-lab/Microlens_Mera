%But : enregistrer les données pour chaque champ d'OGLE IV avec les comparaison avec es expériences
%A faire tourner dans la console après avoir fait tourner main_microlens pour enregistrer si nécéssaire

%chemin d'enregistremen
path_field = strcat('../graph/OGLEIV/', field);

%Création dossier
mkdir(path_field);

%Enregistrement figure efficacité
saveas(1, strcat(path_field, '/efficacité_', field), 'epsc')

%Enregistrement hist normalisé
saveas(18, strcat(path_field, '/hist_normalise_', field), 'epsc')

%Enregistrement hist à l'échelle
% saveas(17, strcat(path_field, '/hist_scale_', field), 'epsc')
% print(17, '-depsc', '-r600', strcat(path_field, '/hist_scale_', field))
% saveas(17, strcat(path_field, '/hist_', field, '.png'))
print(17, '-dpng','-r600', strcat(path_field, '/hist_', field));

%Enregistrement des données
file = fopen(strcat(path_field, '/data_', field),'w');

fprintf(file, strcat('Donnée du champ ', field, ' \n'));
fprintf(file, '(long, lat) = (%5.4f,%5.4f) \n',  table6.glon(table6.field == field), table6.glat(table6.field == field));
fprintf(file, 'N_events = %5.1f \n', table7.N_events(table7.field == field));
fprintf(file, 'N_star = %5.1f \n', table7.N_stars(table7.field == field));
fprintf(file, '\n');
fprintf(file, 'On enregistre des données importantes pour pouvoir les réustiliser ensuite \n');
fprintf(file, 'données expérimentales de OGLE : \n');
fprintf(file, 'gamma_ogle =  %12.8f \n', table7.gam(find(table7.field == field)));
fprintf(file, 'tau_ogle =  %12.8f 10^-6 \n', table7.t_E_mean(find(table7.field == field)));
fprintf(file, '\n');
fprintf(file, 'Données calculées à partid du modèle : \n');
fprintf(file, 'Données sans efficacité expérimentale intégrées par le MC : ');
fprintf(file, 'gamma : %12.8f \n ', gam);
fprintf(file, 'tau : %12.8f 10^-6 \n', tau*1e6);
fprintf(file, '\n');
fprintf(file, 'Avec efficacité expérimentale : ');
fprintf(file, 'gamma avec  eff = %12.8f \n', gamobs);
fprintf(file, 'tau avec eff : %12.8f 10^-6 \n', tauobs*1e6);
fprintf(file, 'gamma avec eff et blending: %12.8f \n', gamobsb);
fprintf(file, 'tau avec eff et blending: %12.8f 10^-6\n', tauobsb*1e6);


