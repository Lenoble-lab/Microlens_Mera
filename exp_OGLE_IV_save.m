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

fprintf(strcat('Donnée du champ \n', field)
fprintf('On enregistre des données importantes pour pouvoir les réustiliser ensuite \n')
fprintf('données expérimentales de OGLE : \n')
fprintf('gamma_ogle =  %12.8f \n', table7.gam(find(table7.field == field)))
fprintf('tau_ogle =  %12.8f \n', table7.t_E_mean(find(table7.field == field)))
fprintf('\n')
fprintf('Données calculées à partid du modèle : \n')
fprintf('Données sans efficacité expérimentale intégrées par le MC : ')
fprintf('gamma : %12.8f \n ', gam)
fprintf('tau : %12.8f \n', tau)
fprintf('\n')
fprintf('Avec efficacité expérimentale : ')
fprintf('gamma avec  eff = %12.8f \n', gamobs)
fprintf('tau avec eff : %12.8f \n', tauobs)
fprintf('gamma avec eff et blending: %12.8f \n', gamobsb)
fprintf('tau avec eff et blending: %12.8f \n', tauobsb)
