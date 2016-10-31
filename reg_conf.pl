use strict;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

my @aliases;


new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
-host =>'localhost',
-user =>'mysql',
-port =>3306,
-pass =>'***',
-species =>'Multi',
-group =>'compara',
-db_version =>'85',
-dbname =>'ensembl_compara_85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'macaca_mulatta', '-group'   => 'core', '-dbname' => 'macaca_mulatta_core_85_10', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'echinops_telfairi', '-group'   => 'core', '-dbname' => 'echinops_telfairi_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'tupaia_belangeri', '-group'   => 'core', '-dbname' => 'tupaia_belangeri_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'erinaceus_europaeus', '-group'   => 'core', '-dbname' => 'erinaceus_europaeus_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'sorex_araneus', '-group'   => 'core', '-dbname' => 'sorex_araneus_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'microcebus_murinus', '-group'   => 'core', '-dbname' => 'microcebus_murinus_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'pongo_abelii', '-group'   => 'core', '-dbname' => 'pongo_abelii_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'equus_caballus', '-group'   => 'core', '-dbname' => 'equus_caballus_core_85_2', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'ochotona_princeps', '-group'   => 'core', '-dbname' => 'ochotona_princeps_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'cavia_porcellus', '-group'   => 'core', '-dbname' => 'cavia_porcellus_core_85_3', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'choloepus_hoffmanni', '-group'   => 'core', '-dbname' => 'choloepus_hoffmanni_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'procavia_capensis', '-group'   => 'core', '-dbname' => 'procavia_capensis_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'tursiops_truncatus', '-group'   => 'core', '-dbname' => 'tursiops_truncatus_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'tarsius_syrichta', '-group'   => 'core', '-dbname' => 'tarsius_syrichta_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'dipodomys_ordii', '-group'   => 'core', '-dbname' => 'dipodomys_ordii_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'vicugna_pacos', '-group'   => 'core', '-dbname' => 'vicugna_pacos_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'pteropus_vampyrus', '-group'   => 'core', '-dbname' => 'pteropus_vampyrus_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'loxodonta_africana', '-group'   => 'core', '-dbname' => 'loxodonta_africana_core_85_3', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'oryctolagus_cuniculus', '-group'   => 'core', '-dbname' => 'oryctolagus_cuniculus_core_85_2', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'ailuropoda_melanoleuca', '-group'   => 'core', '-dbname' => 'ailuropoda_melanoleuca_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'nomascus_leucogenys', '-group'   => 'core', '-dbname' => 'nomascus_leucogenys_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'callithrix_jacchus', '-group'   => 'core', '-dbname' => 'callithrix_jacchus_core_85_321', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'myotis_lucifugus', '-group'   => 'core', '-dbname' => 'myotis_lucifugus_core_85_2', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'bos_taurus', '-group'   => 'core', '-dbname' => 'bos_taurus_core_85_31', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'gorilla_gorilla', '-group'   => 'core', '-dbname' => 'gorilla_gorilla_core_85_31', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'otolemur_garnettii', '-group'   => 'core', '-dbname' => 'otolemur_garnettii_core_85_3', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'pan_troglodytes', '-group'   => 'core', '-dbname' => 'pan_troglodytes_core_85_214', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'ictidomys_tridecemlineatus', '-group'   => 'core', '-dbname' => 'ictidomys_tridecemlineatus_core_85_2', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'sus_scrofa', '-group'   => 'core', '-dbname' => 'sus_scrofa_core_85_102', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'mus_musculus', '-group'   => 'core', '-dbname' => 'mus_musculus_core_85_38', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'canis_familiaris', '-group'   => 'core', '-dbname' => 'canis_familiaris_core_85_31', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'mustela_putorius_furo', '-group'   => 'core', '-dbname' => 'mustela_putorius_furo_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'felis_catus', '-group'   => 'core', '-dbname' => 'felis_catus_core_85_62', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'ovis_aries', '-group'   => 'core', '-dbname' => 'ovis_aries_core_85_31', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'dasypus_novemcinctus', '-group'   => 'core', '-dbname' => 'dasypus_novemcinctus_core_85_3', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'homo_sapiens', '-group'   => 'core', '-dbname' => 'homo_sapiens_core_85_38', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'papio_anubis', '-group'   => 'core', '-dbname' => 'papio_anubis_core_85_2', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'chlorocebus_sabaeus', '-group'   => 'core', '-dbname' => 'chlorocebus_sabaeus_core_85_1', '-db_version'=>'85');

new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'=>'localhost','-user'=> 'mysql', '-port'=> '3306', '-pass'=> '***', '-species' => 'rattus_norvegicus', '-group'   => 'core', '-dbname' => 'rattus_norvegicus_core_85_6', '-db_version'=>'85');
1;