#!/bin/bash
#$ -cwd
#$ -R y
#$ -e wig2bed.err
#$ -o wig2bed.out
#$ -q h14.q
#$ -pe ompi511h14 2
#$ -V                    #export environment var
#$ -N phylop             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

date

export PATH=$PATH:~/programari/cactus-bin-v2.6.8/bin
export PYTHONPATH=~/programari/cactus-bin-v2.6.8/lib:\$PYTHONPATH
source activate py390
source /users-d3/jferrer/programari/cactus-bin-v2.6.8/venv-cactus-v2.6.8/bin/activate

# 1. Get list of those SNPs and also write a file with the SNP and the reference allele

zcat HIGH.vcf.gz | grep -v "#" | cut -f 1,2 | awk -F '\t' '{print $1, $2 - 1, $2}' | sed "s/ /\t/g" > CONS.SNPs.bed
zcat HIGH.vcf.gz | grep -v "#" | cut -f 1,2,4 > CONS.SNPs.csv

# 2. For each SNP, convert to maf and count the number of each base pair; write down the major alelel

touch CONS.SNPs.refallele.csv

while read -r line
do

echo $line | sed "s/ /\t/g" > prova.bed

j=$(cut -f1 prova.bed)
i=$(cut -f2 prova.bed)

halSnps --length 1 --refSequence $j --start $i --tsv prova.tsv ~/gizquierdo/TFM/Cactus/Cactus/363-avian-2020.wPmau.hal Pmau Podiceps_cristatus,Podilymbus_podiceps,Phoenicopterus_ruber,Mesitornis_unicolor,Syrrhaptes_paradoxus,Pterocles_gutturalis,Pterocles_burchelli,Columbina_picui,Patagioenas_fasciata,Columba_livia,Alopecoenas_beccarii,Caloenas_nicobarica,Steatornis_caripensis,Nyctibius_grandis,Nyctibius_bracteatus,Podargus_strigoides,Aegotheles_bennettii,Chaetura_pelagica,Hemiprocne_comata,Oreotrochilus_melanogaster,Calypte_anna,Chordeiles_acutipennis,Antrostomus_carolinensis,Nyctiprogne_leucopyga,Tauraco_erythrolophus,Corythaixoides_concolor,Corythaeola_cristata,Lophotis_ruficrista,Chlamydotis_macqueenii,Ardeotis_kori,Crotophaga_sulcirostris,Geococcyx_californianus,Centropus_bengalensis,Centropus_unirufus,Ceuthmochares_aereus,Piaya_cayana,Cuculus_canorus,Leptosomus_discolor,Rhinopomastus_cyanomelas,Upupa_epops,Buceros_rhinoceros,Bucorvus_abyssinicus,Merops_nubicus,Eurystomus_gularis,Brachypteracias_leptosomus,Todus_mexicanus,Ceyx_cyanopectus,Chloroceryle_aenea,Halcyon_senegalensis,Baryphthengus_martii,Bucco_capensis,Galbula_dea,Picoides_pubescens,Indicator_maculatus,Psilopogon_haemacephalus,Eubucco_bourcierii,Ramphastos_sulfuratus,Semnornis_frantzii,Tricholaema_leucomelas,Trogon_melanurus,Apaloderma_vittatum,Urocolius_indicus,Colius_striatus,Tyto_alba,Strix_occidentalis,Ciccaba_nigrolineata,Glaucidium_brasilianum,Cathartes_aura,Sagittarius_serpentarius,Pandion_haliaetus,Circaetus_pectoralis,Haliaeetus_albicilla,Haliaeetus_leucocephalus,Spizaetus_tyrannus,Aquila_chrysaetos,Cariama_cristata,Chunga_burmeisteri,Nestor_notabilis,Melopsittacus_undulatus,Agapornis_roseicollis,Amazona_guildingii,Probosciger_aterrimus,Eolophus_roseicapillus,Acanthisitta_chloris,Atrichornis_clamosus,Menura_novaehollandiae,Ptilonorhynchus_violaceus,Climacteris_rufus,Malurus_elegans,Dasyornis_broadbenti,Origma_solitaria,Pardalotus_punctatus,Grantiella_picta,Pomatostomus_ruficeps,Orthonyx_spaldingii,Ptilorrhoa_leucosticta,Mohoua_ochrocephala,Daphoenositta_chrysoptera,Eulacestoma_nigropectus,Oreocharis_arfaki,Pteruthius_melanotis,Vireo_altiloquus,Erpornis_zantholeuca,Oriolus_oriolus,Pachycephala_philippinensis,Falcunculus_frontatus,Aleadryas_rufinucha,Chaetorhynchus_papuensis,Rhipidura_dahli,Dicrurus_megarhynchus,Struthidea_cinerea,Aphelocoma_coerulescens,Corvus_moneduloides,Corvus_brachyrhynchos,Corvus_cornix,Lanius_ludovicianus,Ifrita_kowaldi,Paradisaea_raggiana,Myiagra_hebetior,Machaerirhynchus_nigripectus,Rhagologus_leucostigma,Dryoscopus_gambensis,Dyaphorophyia_castanea,Mystacornis_crossleyi,Gymnorhina_tibicen,Edolisoma_coerulescens,Cnemophilus_loriae,Notiomystis_cincta,Callaeas_wilsoni,Anthoscopus_minutus,Poecile_atricapillus,Pseudopodoces_humilis,Parus_major,Panurus_biarmicus,Alaudala_cheleensis,Sylvietta_virens,Cisticola_juncidis,Hippolais_icterina,Acrocephalus_arundinaceus,Oxylabes_madagascariensis,Donacobius_atricapilla,Locustella_ochotensis,Hirundo_rustica,Phylloscopus_trochilus,Rhadina_sibilatrix,Hylia_prasina,Aegithalos_caudatus,Cettia_cetti,Horornis_vulcanius,Erythrocercus_mccallii,Sinosuthora_webbiana,Sylvia_borin,Sylvia_atricapilla,Pomatorhinus_ruficollis,Illadopsis_cleaveri,Leiothrix_lutea,Zosterops_hypoxanthus,Zosterops_lateralis,Sterrhoptilus_dennistouni,Pycnonotus_jocosus,Brachypodius_atriceps,Nicator_chloris,Tichodroma_muraria,Thryothorus_ludovicianus,Polioptila_caerulea,Certhia_familiaris,Certhia_brachydactyla,Sitta_europaea,Elachura_formosa,Cinclus_mexicanus,Buphagus_erythrorhynchus,Toxostoma_redivivum,Leucopsar_rothschildi,Sturnus_vulgaris,Rhabdornis_inornatus,Catharus_fuscescens,Ficedula_albicollis,Saxicola_maurus,Oenanthe_oenanthe,Erithacus_rubecula,Copsychus_sechellarum,Cercotrichas_coryphaeus,Regulus_satrapa,Phainopepla_nitens,Bombycilla_garrulus,Promerops_cafer,Dicaeum_eximium,Leptocoma_aspasia,Chloropsis_cyanopogon,Chloropsis_hardwickii,Peucedramus_taeniatus,Urocynchramus_pylzowi,Ploceus_nigricollis,Vidua_chalybeata,Vidua_macroura,Lonchura_striata,Taeniopygia_guttata,Prunella_fulvescens,Prunella_himalayana,Passer_domesticus,Hypocryptadius_cinnamomeus,Motacilla_alba,Hemignathus_wilsoni,Chlorodrepanis_virens,Serinus_canaria,Loxia_curvirostra,Loxia_leucoptera,Rhodinocichla_rosea,Calcarius_ornatus,Emberiza_fucata,Nesospiza_acunhae,Geospiza_fortis,Sporophila_hypoxantha,Pheucticus_melanocephalus,Cardinalis_cardinalis,Passerina_amoena,Quiscalus_mexicanus,Molothrus_ater,Agelaius_phoeniceus,Setophaga_coronata,Setophaga_kirtlandii,Zonotrichia_albicollis,Junco_hyemalis,Melospiza_melodia,Spizella_passerina,Irena_cyanogastra,Picathartes_gymnocephalus,Chaetops_frenatus,Drymodes_brunneopygia,Melanocharis_versteri,Cephalopterus_ornatus,Pachyramphus_minor,Onychorhynchus_coronatus,Oxyruncus_cristatus,Piprites_chloris,Neopipo_cinnamomea,Tachuris_rubrigastra,Tyrannus_savana,Mionectes_macconnelli,Lepidothrix_coronata,Manacus_manacus,Rhegmatorhina_hoffmannsi,Sakesphorus_luctuosus,Grallaria_varia,Scytalopus_superciliaris,Formicarius_rufipectus,Sclerurus_mexicanus,Furnarius_figulus,Campylorhamphus_procurvoides,Xiphorhynchus_elegans,Serilophus_lunatus,Neodrepanis_coruscans,Pitta_sordida,Sapayoa_aenigma,Calyptomena_viridis,Smithornis_capensis,Herpetotheres_cachinnans,Falco_cherrug,Falco_peregrinus,Fregata_magnificens,Sula_dactylatra,Anhinga_rufa,Anhinga_anhinga,Phalacrocorax_pelagicus,Phalacrocorax_carbo,Phalacrocorax_harrisi,Phalacrocorax_auritus,Phalacrocorax_brasilianus,Egretta_garzetta,Cochlearius_cochlearius,Scopus_umbretta,Pelecanus_crispus,Balaeniceps_rex,Nipponia_nippon,Mesembrinibis_cayennensis,Ciconia_maguari,Aptenodytes_forsteri,Pygoscelis_adeliae,Thalassarche_chlororhynchos,Oceanites_oceanicus,Fregetta_grallaria,Fulmarus_glacialis,Pelecanoides_urinatrix,Calonectris_borealis,Hydrobates_tethys,Gavia_stellata,Phaethon_lepturus,Eurypyga_helias,Rhynochetos_jubatus,Opisthocomus_hoazin,Charadrius_vociferus,Charadrius_alexandrinus,Ibidorhyncha_struthersii,Himantopus_himantopus,Burhinus_bistriatus,Pluvianellus_socialis,Chionis_minor,Pedionomus_torquatus,Thinocorus_orbignyianus,Jacana_jacana,Nycticryphes_semicollaris,Rostratula_benghalensis,Limosa_lapponica,Calidris_pugnax,Arenaria_interpres,Scolopax_mira,Turnix_velox,Dromas_ardeola,Rhinoptilus_africanus,Glareola_pratincola,Rynchops_niger,Phaetusa_simplex,Chroicocephalus_maculipennis,Larus_smithsonianus,Rissa_tridactyla,Stercorarius_parasiticus,Cepphus_grylle,Alca_torda,Uria_lomvia,Uria_aalge,Psophia_crepitans,Grus_americana,Balearica_regulorum,Aramus_guarauna,Heliornis_fulica,Zapornia_atra,Atlantisia_rogersi,Chauna_torquata,Anser_cygnoid,Cairina_moschata,Asarcornis_scutulata,Anas_zonorhyncha,Anas_platyrhynchos,Anseranas_semipalmata,Alectura_lathami,Penelope_pileata,Numida_meleagris,Odontophorus_gujanensis,Callipepla_squamata,Colinus_virginianus,Gallus_gallus,Coturnix_japonica,Phasianus_colchicus,Meleagris_gallopavo,Tympanuchus_cupido,Dromaius_novaehollandiae,Casuarius_casuarius,Apteryx_australis,Apteryx_rowi,Apteryx_haastii,Apteryx_owenii,Nothocercus_nigrocapillus,Nothocercus_julius,Tinamus_guttatus,Crypturellus_soui,Crypturellus_cinnamomeus,Crypturellus_undulatus,Nothoprocta_perdicaria,Nothoprocta_ornata,Nothoprocta_pentlandii,Eudromia_elegans,Rhea_americana,Rhea_pennata,Struthio_camelus

countG=$(less prova.tsv | tail -n -1 | grep -o -n -e "G" -e "g" | wc -l)
countA=$(less prova.tsv | tail -n -1 | grep -o -n -e "A" -e "A" | wc -l)
countC=$(less prova.tsv | tail -n -1 | grep -o -n -e "C" -e "c" | wc -l)
countT=$(less prova.tsv | tail -n -1 | grep -o -n -e "T" -e "t" | wc -l)

echo -e "G\t${countG}" > bases.csv
echo -e "A\t${countA}" >> bases.csv
echo -e "C\t${countC}" >> bases.csv
echo -e "T\t${countT}" >> bases.csv

sort -k2 -n -r bases.csv | head -n 1 | cut -f 1 >> CONS.SNPs.refallele.csv

rm bases.csv
rm prova.bed

done < CONS.SNPs.bed

paste CONS.SNPs.csv CONS.SNPs.refallele.csv > prova
mv prova CONS.SNPs.csv

# 3. Subset VCFs with REFdel and ALTdel poisitions

cat CONS.SNPs.csv | awk '$3 == $4' | cut -f 1,2 > Puffinus.SNP.ALTdel.bed
cat CONS.SNPs.csv | awk '$3 != $4' | cut -f 1,2 > Puffinus.SNP.REFdel.bed

vcftools --gzvcf HIGH.vcf.gz --positions Puffinus.SNP.ALTdel.bed --recode-INFO-all --recode --out Puffinus.SNP.ALTdel
mv Puffinus.SNP.ALTdel.recode.vcf Puffinus.SNP.ALTdel.vcf
bgzip Puffinus.SNP.ALTdel.vcf
tabix Puffinus.SNP.ALTdel.vcf.gz

vcftools --gzvcf HIGH.vcf.gz --positions Puffinus.SNP.REFdel.bed --recode-INFO-all --recode --out Puffinus.SNP.REFdel
mv Puffinus.SNP.REFdel.recode.vcf Puffinus.SNP.REFdel.vcf
bgzip Puffinus.SNP.REFdel.vcf
tabix Puffinus.SNP.REFdel.vcf.gz

# 4. Repolarize the REFdel dataset and join both in a single file

gunzip Puffinus.SNP.REFdel.vcf.gz

cat Puffinus.SNP.REFdel.vcf | sed 's/0\/0/0|0/g' > prova
mv prova Puffinus.SNP.REFdel.vcf
cat Puffinus.SNP.REFdel.vcf | sed 's/1\/1/1|1/g' > prova
mv prova Puffinus.SNP.REFdel.vcf
cat Puffinus.SNP.REFdel.vcf | sed 's/1|1/0\/0/g' > prova
mv prova Puffinus.SNP.REFdel.vcf
cat Puffinus.SNP.REFdel.vcf | sed 's/0|0/1\/1/g' > prova
mv prova Puffinus.SNP.REFdel.vcf

bgzip Puffinus.SNP.REFdel.vcf

bcftools concat -Oz -o Puffinus.SNP.polarized.vcf.gz -a --threads 12 Puffinus.SNP.REFdel.vcf.gz Puffinus.SNP.ALTdel.vcf.gz
tabix Puffinus.SNP.polarized.vcf.gz

###### THEN FOLLOW THE GUIDELINES IN ../Cactus&phyloP/7_calculate_load.sh
