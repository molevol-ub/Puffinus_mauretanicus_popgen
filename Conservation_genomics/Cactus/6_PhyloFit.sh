#!/bin/bash
#$ -cwd
#$ -R y
#$ -e phyloFit.err
#$ -o phyloFit.out
#$ -q h13.q
#$ -pe ompi255h13 2
#$ -V                    #export environment var
#$ -N phyloFit             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

date

export PATH=$PATH:~/programari/cactus-bin-v2.5.1/bin
source activate py39
cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus

# 1. Convert hal to maf



# 2. Run phyloFit
# Mostly following Patrick F. Sullivan et al. (2023) - Science on mammals  

phyloFit --tree "((((((Podiceps_cristatus:0.0392987,Podilymbus_podiceps:0.0235222)birdAnc0:0.0429797,Phoenicopterus_ruber:0.0405362)birdAnc1:0.0167344,((Mesitornis_unicolor:0.10254,((Syrrhaptes_paradoxus:0.0425251,Pterocles_gutturalis:0.0262198)birdAnc2:0.00460539,Pterocles_burchelli:0.0284758)birdAnc3:0.0527532)birdAnc4:0.00537808,(Columbina_picui:0.0453656,((Patagioenas_fasciata:0.0257823,Columba_livia:0.0273616)birdAnc5:0.0192441,(Alopecoenas_beccarii:0.0393388,Caloenas_nicobarica:0.0323221)birdAnc6:0.00386343)birdAnc7:0.00308433)birdAnc8:0.0741755)birdAnc9:0.00549477)birdAnc10:0.000962904,(((((Steatornis_caripensis:0.0666685,(Nyctibius_grandis:0.0253729,Nyctibius_bracteatus:0.02909)birdAnc11:0.0344433)birdAnc12:0.00504245,(Podargus_strigoides:0.0947154,(Aegotheles_bennettii:0.0798626,((Chaetura_pelagica:0.0634338,Hemiprocne_comata:0.0399419)birdAnc13:0.0294833,(Oreotrochilus_melanogaster:0.0256801,Calypte_anna:0.0221997)birdAnc14:0.0929121)birdAnc15:0.0176493)birdAnc16:0.0105239)birdAnc17:0.00613956)birdAnc18:0.00188528,((Chordeiles_acutipennis:0.0340503,Antrostomus_carolinensis:0.0339401)birdAnc19:0.00414229,Nyctiprogne_leucopyga:0.0444645)birdAnc20:0.0488615)birdAnc21:0.00667159,((Tauraco_erythrolophus:0.0541986,(Corythaixoides_concolor:0.0287367,Corythaeola_cristata:0.0392096)birdAnc22:0.0144931)birdAnc23:0.0382012,((Lophotis_ruficrista:0.0277471,(Chlamydotis_macqueenii:0.0333016,Ardeotis_kori:0.0218211)birdAnc24:0.00511901)birdAnc25:0.0506743,((Crotophaga_sulcirostris:0.0544226,Geococcyx_californianus:0.0675258)birdAnc26:0.0146249,((Centropus_bengalensis:0.0205289,Centropus_unirufus:0.0136295)birdAnc27:0.0910811,((Ceuthmochares_aereus:0.0394053,Piaya_cayana:0.059911)birdAnc28:0.0189657,Cuculus_canorus:0.0473037)birdAnc29:0.0226993)birdAnc30:0.0123658)birdAnc31:0.0403424)birdAnc32:0.00521318)birdAnc33:0.00317325)birdAnc34:0.00449453,((((((Leptosomus_discolor:0.0842217,((((Rhinopomastus_cyanomelas:0.0998847,Upupa_epops:0.0951532)birdAnc35:0.0885654,(Buceros_rhinoceros:0.057707,Bucorvus_abyssinicus:0.030074)birdAnc36:0.0513299)birdAnc37:0.0141845,(((Merops_nubicus:0.164824,(Eurystomus_gularis:0.0554788,Brachypteracias_leptosomus:0.0899502)birdAnc38:0.0204202)birdAnc39:0.0069423,(Todus_mexicanus:0.13241,((Ceyx_cyanopectus:0.179304,(Chloroceryle_aenea:0.0689874,Halcyon_senegalensis:0.0457873)birdAnc40:0.00640787)birdAnc41:0.0181018,Baryphthengus_martii:0.0946033)birdAnc42:0.0112638)birdAnc43:0.00865139)birdAnc44:0.0148702,((Bucco_capensis:0.0960093,Galbula_dea:0.0836842)birdAnc45:0.0269525,((Picoides_pubescens:0.0577799,Indicator_maculatus:0.0421248)birdAnc46:0.0267434,(Psilopogon_haemacephalus:0.0754123,((Eubucco_bourcierii:0.0308309,(Ramphastos_sulfuratus:0.0238005,Semnornis_frantzii:0.0533906)birdAnc47:0.0082454)birdAnc48:0.0287688,Tricholaema_leucomelas:0.0557189)birdAnc49:0.00460952)birdAnc50:0.0317959)birdAnc51:0.0795069)birdAnc52:0.0210906)birdAnc53:0.00780139)birdAnc54:0.00602991,(Trogon_melanurus:0.0358249,Apaloderma_vittatum:0.0458747)birdAnc55:0.0793133)birdAnc56:0.00900188)birdAnc57:0.00646619,(Urocolius_indicus:0.046402,Colius_striatus:0.0467654)birdAnc58:0.101228)birdAnc59:0.00926588,((Tyto_alba:0.0713578,((Strix_occidentalis:0.00447592,Ciccaba_nigrolineata:0.00601834)birdAnc60:0.0158527,Glaucidium_brasilianum:0.0402936)birdAnc61:0.0326469)birdAnc62:0.0129786,(Cathartes_aura:0.0504201,(Sagittarius_serpentarius:0.0523745,(Pandion_haliaetus:0.0522228,(Circaetus_pectoralis:0.0134599,((Haliaeetus_albicilla:0.0011328,Haliaeetus_leucocephalus:0.00824849)birdAnc63:0.0163168,(Spizaetus_tyrannus:0.00861222,Aquila_chrysaetos:0.0190887)birdAnc64:0.00681741)birdAnc65:0.00268443)birdAnc66:0.019219)birdAnc67:0.0109663)birdAnc68:0.0080475)birdAnc69:0.0198074)birdAnc70:0.00114941)birdAnc71:0.00109145,((Cariama_cristata:0.0154642,Chunga_burmeisteri:0.0104022)birdAnc72:0.0512415,(((Nestor_notabilis:0.0476696,(((Melopsittacus_undulatus:0.0385032,Agapornis_roseicollis:0.0454335)birdAnc73:0.00900727,Amazona_guildingii:0.0403314)birdAnc74:0.00681355,(Probosciger_aterrimus:0.0258601,Eolophus_roseicapillus:0.0185469)birdAnc75:0.0166148)birdAnc76:0.0154575)birdAnc77:0.0510807,(Acanthisitta_chloris:0.102918,(((Atrichornis_clamosus:0.0421806,Menura_novaehollandiae:0.034916)birdAnc78:0.0192208,((Ptilonorhynchus_violaceus:0.0466603,Climacteris_rufus:0.0764537)birdAnc79:0.0125135,((Malurus_elegans:0.0719273,(Dasyornis_broadbenti:0.0474528,((Origma_solitaria:0.0431673,Pardalotus_punctatus:0.0396993)birdAnc80:0.00671108,Grantiella_picta:0.0525609)birdAnc81:0.00110922)birdAnc82:0.00396686)birdAnc83:0.0119031,((Pomatostomus_ruficeps:0.0754387,Orthonyx_spaldingii:0.0454468)birdAnc84:0.00762959,((Ptilorrhoa_leucosticta:0.0385926,(((Mohoua_ochrocephala:0.0451064,Daphoenositta_chrysoptera:0.0487468)birdAnc85:0.00785046,((Eulacestoma_nigropectus:0.0319099,(((Oreocharis_arfaki:0.0343827,(Pteruthius_melanotis:0.0337533,(Vireo_altiloquus:0.034824,Erpornis_zantholeuca:0.0353804)birdAnc86:0.00990282)birdAnc87:0.00809337)birdAnc88:0.00102241,(Oriolus_oriolus:0.0318345,Pachycephala_philippinensis:0.0349829)birdAnc89:0.00618331)birdAnc90:0.0013154,(Falcunculus_frontatus:0.0273883,Aleadryas_rufinucha:0.0331317)birdAnc91:0.00970177)birdAnc92:0.00127741)birdAnc93:0.00121279,(((Chaetorhynchus_papuensis:0.0270311,Rhipidura_dahli:0.0350931)birdAnc94:0.0103699,(Dicrurus_megarhynchus:0.0272301,((Struthidea_cinerea:0.0274171,((Aphelocoma_coerulescens:0.0140808,(Corvus_moneduloides:0.01015,(Corvus_brachyrhynchos:0.003901,Corvus_cornix:0.00224804)birdAnc95:0.0026649)birdAnc96:0.00982717)birdAnc97:0.0123772,Lanius_ludovicianus:0.0372667)birdAnc98:0.00442354)birdAnc99:0.00120252,((Ifrita_kowaldi:0.027878,Paradisaea_raggiana:0.021846)birdAnc100:0.00467297,Myiagra_hebetior:0.0391168)birdAnc101:0.00101262)birdAnc102:0.00131007)birdAnc103:0.00129648)birdAnc104:0.00794971,(Machaerirhynchus_nigripectus:0.0316523,((Rhagologus_leucostigma:0.0267596,(Dryoscopus_gambensis:0.0357394,(Dyaphorophyia_castanea:0.0269874,Mystacornis_crossleyi:0.0468135)birdAnc105:0.00327938)birdAnc106:0.0023785)birdAnc107:0.0051592,Gymnorhina_tibicen:0.0226897)birdAnc108:0.00129249)birdAnc109:0.00626371)birdAnc110:0.00118305)birdAnc111:0.000329723)birdAnc112:0.000910678,Edolisoma_coerulescens:0.0329734)birdAnc113:0.00257603)birdAnc114:0.00570099,(Cnemophilus_loriae:0.0400785,(((Notiomystis_cincta:0.0732189,Callaeas_wilsoni:0.036756)birdAnc115:0.00817619,((((Anthoscopus_minutus:0.0502255,(Poecile_atricapillus:0.0168745,(Pseudopodoces_humilis:0.0109522,Parus_major:0.0128199)birdAnc116:0.00505373)birdAnc117:0.0253478)birdAnc118:0.00860147,((Panurus_biarmicus:0.0555984,Alaudala_cheleensis:0.0785207)birdAnc119:0.00980353,((Sylvietta_virens:0.0416221,(Cisticola_juncidis:0.0808908,(((Hippolais_icterina:0.0227433,Acrocephalus_arundinaceus:0.0251951)birdAnc120:0.0248349,((Oxylabes_madagascariensis:0.0640696,Donacobius_atricapilla:0.046076)birdAnc121:0.00860278,Locustella_ochotensis:0.0451303)birdAnc122:0.00647972)birdAnc123:0.00248291,(Hirundo_rustica:0.0532979,(((Phylloscopus_trochilus:0.0300936,Rhadina_sibilatrix:0.0266193)birdAnc124:0.0210544,(Hylia_prasina:0.0432648,(Aegithalos_caudatus:0.0561489,((Cettia_cetti:0.0358456,Horornis_vulcanius:0.0364962)birdAnc125:0.0157272,Erythrocercus_mccallii:0.0466275)birdAnc126:0.00303645)birdAnc127:0.00131804)birdAnc128:0.00286041)birdAnc129:0.00606578,(((Sinosuthora_webbiana:0.0408785,(Sylvia_borin:0.0149836,Sylvia_atricapilla:0.0236626)birdAnc130:0.0236187)birdAnc131:0.00572285,((Pomatorhinus_ruficollis:0.0325245,(Illadopsis_cleaveri:0.0284669,Leiothrix_lutea:0.0461797)birdAnc132:0.00591997)birdAnc133:0.00319567,((Zosterops_hypoxanthus:0.00737955,Zosterops_lateralis:0.00600487)birdAnc134:0.0131015,Sterrhoptilus_dennistouni:0.0174757)birdAnc135:0.0350534)birdAnc136:0.00385827)birdAnc137:0.0121306,(Pycnonotus_jocosus:0.0274686,Brachypodius_atriceps:0.0296084)birdAnc138:0.0281225)birdAnc139:0.00521403)birdAnc140:0.00259651)birdAnc141:0.00104268)birdAnc142:0.00105669)birdAnc143:0.00598171)birdAnc144:0.00461601,Nicator_chloris:0.0501384)birdAnc145:0.00115739)birdAnc146:0.004682)birdAnc147:0.00160681,((((Tichodroma_muraria:0.0601018,(((Thryothorus_ludovicianus:0.0453808,Polioptila_caerulea:0.039363)birdAnc148:0.0190269,(Certhia_familiaris:0.0129074,Certhia_brachydactyla:0.00587503)birdAnc149:0.0408078)birdAnc150:0.00286194,Sitta_europaea:0.0563427)birdAnc151:0.0018102)birdAnc152:0.0056281,(Elachura_formosa:0.0824017,(Cinclus_mexicanus:0.0515478,((Buphagus_erythrorhynchus:0.0416008,(Toxostoma_redivivum:0.0306843,((Leucopsar_rothschildi:0.00966742,Sturnus_vulgaris:0.00943466)birdAnc153:0.0173548,Rhabdornis_inornatus:0.0251144)birdAnc154:0.00543945)birdAnc155:0.00527706)birdAnc156:0.0135923,(Catharus_fuscescens:0.0506154,(((Ficedula_albicollis:0.0195583,(Saxicola_maurus:0.0150623,Oenanthe_oenanthe:0.0175955)birdAnc157:0.0126836)birdAnc158:0.00467272,Erithacus_rubecula:0.0357394)birdAnc159:0.0144944,(Copsychus_sechellarum:0.0177428,Cercotrichas_coryphaeus:0.019971)birdAnc160:0.0115204)birdAnc161:0.0155683)birdAnc162:0.00615042)birdAnc163:0.00467074)birdAnc164:0.0100852)birdAnc165:0.00306837)birdAnc166:0.00121038,(Regulus_satrapa:0.0724478,(Phainopepla_nitens:0.0339186,Bombycilla_garrulus:0.0414967)birdAnc167:0.0247518)birdAnc168:0.00826091)birdAnc169:0.00202142,(Promerops_cafer:0.038775,((Dicaeum_eximium:0.0438166,Leptocoma_aspasia:0.0283009)birdAnc170:0.012609,((Chloropsis_cyanopogon:0.0256733,Chloropsis_hardwickii:0.0203012)birdAnc171:0.0310394,((Peucedramus_taeniatus:0.0409992,((Urocynchramus_pylzowi:0.0420732,(Ploceus_nigricollis:0.0250551,((Vidua_chalybeata:0.00778432,Vidua_macroura:0.00906249)birdAnc172:0.0189881,(Lonchura_striata:0.0188093,Taeniopygia_guttata:0.0180254)birdAnc173:0.0182355)birdAnc174:0.00632893)birdAnc175:0.00521183)birdAnc176:0.00113078,((Prunella_fulvescens:0.0188495,Prunella_himalayana:0.0120033)birdAnc177:0.0316134,((Passer_domesticus:0.0404815,Hypocryptadius_cinnamomeus:0.0374908)birdAnc178:0.00843521,(Motacilla_alba:0.0389634,(((Hemignathus_wilsoni:0.0117011,Chlorodrepanis_virens:0.00460529)birdAnc179:0.0141782,(Serinus_canaria:0.0137162,(Loxia_curvirostra:0.00532379,Loxia_leucoptera:0.00206904)birdAnc180:0.00761149)birdAnc181:0.010181)birdAnc182:0.0130781,(Rhodinocichla_rosea:0.0218409,((Calcarius_ornatus:0.0297823,Emberiza_fucata:0.0273423)birdAnc183:0.00705278,(((Nesospiza_acunhae:0.0179732,(Geospiza_fortis:0.0198689,Sporophila_hypoxantha:0.0242787)birdAnc184:0.001094)birdAnc185:0.0068089,(Pheucticus_melanocephalus:0.0244591,(Cardinalis_cardinalis:0.0171086,Passerina_amoena:0.0238962)birdAnc186:0.00191497)birdAnc187:0.0054246)birdAnc188:0.00914053,((((Quiscalus_mexicanus:0.0177607,Molothrus_ater:0.00526022)birdAnc189:0.00174671,Agelaius_phoeniceus:0.0113643)birdAnc190:0.0131346,(Setophaga_coronata:0.0116748,Setophaga_kirtlandii:0.00890108)birdAnc191:0.0198616)birdAnc192:0.00334569,(((Zonotrichia_albicollis:0.014764,Junco_hyemalis:0.0126042)birdAnc193:0.00728693,Melospiza_melodia:0.0164274)birdAnc194:0.00676442,Spizella_passerina:0.0570829)birdAnc195:0.016997)birdAnc196:0.006714)birdAnc197:0.00133576)birdAnc198:0.00340184)birdAnc199:0.014975)birdAnc200:0.00952034)birdAnc201:0.0044316)birdAnc202:0.00526676)birdAnc203:0.00441415)birdAnc204:0.00299822)birdAnc205:0.00523039,Irena_cyanogastra:0.0540942)birdAnc206:0.00161269)birdAnc207:0.00311762)birdAnc208:0.0053686)birdAnc209:0.00364265)birdAnc210:0.0022546)birdAnc211:0.0082532,((Picathartes_gymnocephalus:0.0362297,Chaetops_frenatus:0.0473373)birdAnc212:0.0174355,Drymodes_brunneopygia:0.0687406)birdAnc213:0.00649534)birdAnc214:0.00127266)birdAnc215:0.00307299,Melanocharis_versteri:0.042484)birdAnc216:0.00379532)birdAnc217:0.00512779)birdAnc218:0.0120392)birdAnc219:0.00517287)birdAnc220:0.0153146)birdAnc221:0.0148454)birdAnc222:0.0315882,((((Cephalopterus_ornatus:0.0393414,(Pachyramphus_minor:0.0401089,((Onychorhynchus_coronatus:0.0387256,Oxyruncus_cristatus:0.0334472)birdAnc223:0.00261503,(Piprites_chloris:0.0418464,(Neopipo_cinnamomea:0.0360794,(Tachuris_rubrigastra:0.0390523,(Tyrannus_savana:0.0449023,Mionectes_macconnelli:0.0401292)birdAnc224:0.00470013)birdAnc225:0.00133182)birdAnc226:0.00298921)birdAnc227:0.00716689)birdAnc228:0.00327947)birdAnc229:0.00210344)birdAnc230:0.00414044,(Lepidothrix_coronata:0.0128375,Manacus_manacus:0.0143885)birdAnc231:0.0239822)birdAnc232:0.0323095,((Rhegmatorhina_hoffmannsi:0.0283977,Sakesphorus_luctuosus:0.0245601)birdAnc233:0.046322,(Grallaria_varia:0.0610732,(Scytalopus_superciliaris:0.0534954,(Formicarius_rufipectus:0.0689979,(Sclerurus_mexicanus:0.0306808,(Furnarius_figulus:0.0346674,(Campylorhamphus_procurvoides:0.018145,Xiphorhynchus_elegans:0.0265301)birdAnc234:0.0203852)birdAnc235:0.00691259)birdAnc236:0.0180958)birdAnc237:0.00459733)birdAnc238:0.0179597)birdAnc239:0.00579843)birdAnc240:0.0173701)birdAnc241:0.0131603,((Serilophus_lunatus:0.0521748,Neodrepanis_coruscans:0.0839779)birdAnc242:0.0248413,((Pitta_sordida:0.0774412,Sapayoa_aenigma:0.0624119)birdAnc243:0.0068565,(Calyptomena_viridis:0.0494551,Smithornis_capensis:0.048195)birdAnc244:0.0136253)birdAnc245:0.00326511)birdAnc246:0.0284205)birdAnc247:0.0153591)birdAnc248:0.012502)birdAnc249:0.0439386)birdAnc250:0.00913236,(Herpetotheres_cachinnans:0.0290457,(Falco_cherrug:0.003749,Falco_peregrinus:0.00366728)birdAnc251:0.0483553)birdAnc252:0.0281468)birdAnc253:0.0140244)birdAnc254:0.0158096)birdAnc255:0.00624113,((((((Fregata_magnificens:0.0508025,(Sula_dactylatra:0.0538305,((Anhinga_rufa:0.078321,Anhinga_anhinga:0.0235052)birdAnc256:0.0160972,((Phalacrocorax_pelagicus:0.00966167,Phalacrocorax_carbo:0.0114102)birdAnc257:0.00133108,(Phalacrocorax_harrisi:0.00982744,(Phalacrocorax_auritus:0.0125589,Phalacrocorax_brasilianus:0.017429)birdAnc258:0.00297035)birdAnc259:0.0034775)birdAnc260:0.0433022)birdAnc261:0.0064345)birdAnc262:0.0175105)birdAnc263:0.0107208,(((Egretta_garzetta:0.0319996,Cochlearius_cochlearius:0.0256433)birdAnc264:0.0326248,((Scopus_umbretta:0.0857843,Pelecanus_crispus:0.0405858)birdAnc265:0.00116402,Balaeniceps_rex:0.0594473)birdAnc266:0.016726)birdAnc267:0.00349971,(Nipponia_nippon:0.0269239,Mesembrinibis_cayennensis:0.0251115)birdAnc268:0.0198646)birdAnc269:0.0101279)birdAnc270:0.00169584,Ciconia_maguari:0.0401157)birdAnc271:0.00124487,((Aptenodytes_forsteri:0.0115539,Pygoscelis_adeliae:0.016288)birdAnc272:0.024277,(Thalassarche_chlororhynchos:0.0539544,((Oceanites_oceanicus:0.0177612,Fregetta_grallaria:0.0203341)birdAnc273:0.0230702,((Fulmarus_glacialis:0.024219,(Pelecanoides_urinatrix:0.0287616,Calonectris_borealis:0.0270643)birdAnc274:0.00326066)birdAnc275:0.0118373,Hydrobates_tethys:0.0339242)birdAnc276:0.00394151)birdAnc277:0.00280891)birdAnc278:0.00514795)birdAnc279:0.00221182)birdAnc280:0.00135827,Gavia_stellata:0.0572347)birdAnc281:0.00814884,(Phaethon_lepturus:0.0644989,(Eurypyga_helias:0.0740253,Rhynochetos_jubatus:0.0437391)birdAnc282:0.0591763)birdAnc283:0.00639553)birdAnc284:0.00115719)birdAnc285:0.00855485,(Opisthocomus_hoazin:0.092196,(((((Charadrius_vociferus:0.0293659,Charadrius_alexandrinus:0.0256289)birdAnc286:0.0103449,(Ibidorhyncha_struthersii:0.0250383,Himantopus_himantopus:0.0202267)birdAnc287:0.0113403)birdAnc288:0.0152746,(Burhinus_bistriatus:0.0403873,(Pluvianellus_socialis:0.0260733,Chionis_minor:0.0189993)birdAnc289:0.0453714)birdAnc290:0.00968586)birdAnc291:0.00560493,((((Pedionomus_torquatus:0.0571575,Thinocorus_orbignyianus:0.0596294)birdAnc292:0.0137293,(Jacana_jacana:0.0727467,(Nycticryphes_semicollaris:0.030227,Rostratula_benghalensis:0.0439908)birdAnc293:0.0183339)birdAnc294:0.00879479)birdAnc295:0.0156251,(Limosa_lapponica:0.0480545,(Calidris_pugnax:0.0397285,Arenaria_interpres:0.0305997)birdAnc296:0.0112865)birdAnc297:0.01698)birdAnc298:0.0297793,(Turnix_velox:0.194348,((Dromas_ardeola:0.0333739,(Rhinoptilus_africanus:0.0425364,Glareola_pratincola:0.032094)birdAnc299:0.00573402)birdAnc300:0.00864775,(((Rynchops_niger:0.0274225,Phaetusa_simplex:0.0191109)birdAnc301:0.00135514,((Chroicocephalus_maculipennis:0.00448431,Larus_smithsonianus:0.00692431)birdAnc302:0.00775055,Rissa_tridactyla:0.00791143)birdAnc303:0.0151582)birdAnc304:0.0065009,(Stercorarius_parasiticus:0.022584,(Cepphus_grylle:0.0171177,(Alca_torda:0.0105477,(Uria_lomvia:0.0129748,Uria_aalge:0.0133819)birdAnc305:0.00628811)birdAnc306:0.00549562)birdAnc307:0.00934918)birdAnc308:0.00270567)birdAnc309:0.010142)birdAnc310:0.0128998)birdAnc311:0.00618483)birdAnc312:0.0112997)birdAnc313:0.0122781,((Psophia_crepitans:0.0615905,((Grus_americana:0.0221357,Balearica_regulorum:0.0230559)birdAnc314:0.0174302,Aramus_guarauna:0.0329495)birdAnc315:0.0163615)birdAnc316:0.0072558,(Heliornis_fulica:0.0819623,(Zapornia_atra:0.033309,Atlantisia_rogersi:0.0482957)birdAnc317:0.0483444)birdAnc318:0.0355748)birdAnc319:0.0182295)birdAnc320:0.00134879)birdAnc321:0.00113226)birdAnc322:0.0113373)birdAnc323:0.00219801)birdAnc324:0.0363168,((Chauna_torquata:0.05918,((Anser_cygnoid:0.0260992,(Cairina_moschata:0.0163315,(Asarcornis_scutulata:0.0166641,(Anas_zonorhyncha:0.0145983,Anas_platyrhynchos:0.00370956)birdAnc325:0.0101879)birdAnc326:0.00442028)birdAnc327:0.0391372)birdAnc328:0.0622062,Anseranas_semipalmata:0.043174)birdAnc329:0.00628358)birdAnc330:0.00907256,(Alectura_lathami:0.104542,(Penelope_pileata:0.0860068,(Numida_meleagris:0.047111,((Odontophorus_gujanensis:0.0273205,(Callipepla_squamata:0.0144011,Colinus_virginianus:0.0163078)birdAnc331:0.0295101)birdAnc332:0.0438346,((Gallus_gallus:0.0430818,Coturnix_japonica:0.0626533)birdAnc333:0.00650237,(Phasianus_colchicus:0.028144,(Meleagris_gallopavo:0.032322,Tympanuchus_cupido:0.0366513)birdAnc334:0.00419078)birdAnc335:0.0188335)birdAnc336:0.0177949)birdAnc337:0.0046239)birdAnc338:0.0474867)birdAnc339:0.0143208)birdAnc340:0.033988)birdAnc341:0.0360553)birdAnc342:0.0420357,((((Dromaius_novaehollandiae:0.0169367,Casuarius_casuarius:0.0137762)birdAnc343:0.0287306,((Apteryx_australis:0.00625179,Apteryx_rowi:0.00536088)birdAnc344:0.0107418,(Apteryx_haastii:0.0013104,Apteryx_owenii:0.00131367)birdAnc345:0.00428185)birdAnc346:0.0307378)birdAnc347:0.00807181,((((Nothocercus_nigrocapillus:0.0183726,Nothocercus_julius:0.0233871)birdAnc348:0.0326129,(Tinamus_guttatus:0.0466569,(Crypturellus_soui:0.0273277,(Crypturellus_cinnamomeus:0.00965895,Crypturellus_undulatus:0.013486)birdAnc349:0.0142709)birdAnc350:0.0367473)birdAnc351:0.0126899)birdAnc352:0.00679526,((Nothoprocta_perdicaria:0.0208291,(Nothoprocta_ornata:0.0075023,Nothoprocta_pentlandii:0.00824407)birdAnc353:0.0163139)birdAnc354:0.0791579,Eudromia_elegans:0.0798405)birdAnc355:0.00607447)birdAnc356:0.114797,(Rhea_americana:0.0080353,Rhea_pennata:0.00645012)birdAnc357:0.0637885)birdAnc358:0.00609096)birdAnc359:0.0146037,Struthio_camelus:0.0558674)birdAnc360:0.0420357)" \
  --subst-mod REV --EM --out-root 363-avian-2020 363-avian-2020.maf
