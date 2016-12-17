mods = {}
mods['E'] = {        "Phenyl Ester"                             :76,
                     "pyroglutamic acid"                        :-18,
                     "decarboxylation of gamma carboxy"         :-44,
                     "sodium"                                   :22,
                     "Potassium"                                :38,
                     "carboxylation"                            :43,     # vitamin K used to catalyze, mostly in coagulation factors and bone metabolism
                     "Glycerol Ester"                           :74,
                     "O-ADP-ribosylation"                       :541}

mods['C'] = {        "lysinoalanine"                            :-34,
                     "lanthionine"                              :-34,
                     "dehydroalanine"                           :-34,
                     "formylglycine"                            :-18,
                     "sulfenic acid"                            :16,
                     "carbamylation"                            :43,
                     "oxidation to cysteic acid"                :48,
                     "carboxyamidomethyl"                       :57,
                     "carboxymethyl"                            :58,
                     "acetamidomethyl"                          :71,
                     "propionamide adduct (acrylamide)"         :71,
                     "3,5-dichlorination with 37Cl"             :72,
                     "S-(sn-1-Glyceryl)"                        :74,
                     "Beta Mercaptoethanol"                     :76,
                     "pyridylethylation"                        :105,
                     "cysteinylation"                           :119,
                     "farnesylation"                            :204,
                     "palmitoylation"                           :238,
                     "stearoylation"                            :266,
                     "geranylation"                             :272,
                     "S-geranylgeranylation"                    :276,
                     "glutationation"                           :305,
                     "S-(sn-1-Dipalmitoyl-glyceryl"             :541,
                     "S-Phycocyanobilin"                        :587,
                     "S-Heme"                                   :617,
                     "S-(6-Flavin [FAD])"                       :784}

mods['K'] = {        "desmosine"                                :-58,
                     "allysine"                                 :-1,
                     "oxidation to aminoadipic semialdehyde"    :-1,
                     "epsilon amino to imine"                   :12,
                     "syndesine"                                :13,
                     "delta-Hydroxy-allysine"                   :14,
                     "N-epsilon methylation (norleucine)"       :14,
                     "delta-Hydroxy-allysine"                   :15,
                     "oxidation to aminoapidic acid"            :15,
                     "hydroxylation of delta-C"                 :16,
                     "formylation"                              :28,
                     "N epsilon acetylation"                    :42,
                     "N trimethylation"                         :42,
                     "carbamylation"                            :43,
                     "Phosphorylation"                          :80,
                     "phosphate/sulphate adduct"                :98,
                     "pentose addition (Ara, Rib, Xyl)"         :132,
                     "hexose addition (Fru, Gal, Glc, Man)"     :132,
                     "delta-glycosyloxy"                        :177,
                     "lipoic acid"                              :188,
                     "myristoylation"                           :210,
                     "biotinylation"                            :226,
                     "N-pyridoxyl"                              :229,
                     "Pyridoxal phosphate schiff base"          :231}

mods['M'] = {        "decomposed carboxymethylation"            :-48,
                     "homoserine with CNBr"                     :-30,
                     "misincorporation of norleucine"           :-18,
                     "reduction to homocysteine"                :-14,
                     "oxidation to sulphoxide"                  :16,
                     "oxidation to sulphone"                    :32,
                     "Homoseryl lactone with cyanogen bromide"  :83}

mods['R'] = {        "gamma-glutamyl semiladehyde"              :-43,
                     "ornithine"                                :-42,
                     "oxidation of R to E"                      :-27,
                     "citruline"                                :1,
                     "N,N dimethylation"                        :28,
                     "carbamylation"                            :43,
                     "phosphate/sulphate adduct"                :98,
                     "pentose addition (Ara, Rib, Xyl)"         :132,
                     "hexose addition (Fru, Gal, Glc, Man)"     :132,
                     "N-(ADP-ribosyl)"                          :541}

mods['crosslinks'] = {}

mods['crosslinks']['CC'] = {"disulphide bond formation"                :-2,
                            "disulfide reduction"                      :2}

mods['crosslinks']['M0'] = {"N-pyrrolidone carboxyl"                   :-17,
                            "S-carbamoylmethylcysteine cyclization"    :-17}

mods['crosslinks']['YY'] = {"3,3p BiTyr"                               :-2,
                            "IsodiTyr"                                 :-2}

mods['crosslinks']['C0'] = {"aspartate to succinimide"                 :-18}

mods['crosslinks']['CH'] = {"S-2-histidyl xlink to cys"                :-2}
mods['crosslinks']['CY'] = {"S-3-tyr xlink to cys"                     :-2}
mods['crosslinks']['C0'] = {"formaldehyde adduct"                      :12}
mods['crosslinks']['CE'] = {"S-gamma-glutamyl crosslink"               :-18}
mods['crosslinks']['DK'] = {"N-beta-D crosslink"                       :-17}
mods['crosslinks']['EK'] = {"N-alpha-(gamma-glutamyl) crosslink"       :-17}
mods['crosslinks']['ES'] = {"O-gamma-glutamyl crosslink"               :-18}
mods['crosslinks']['HR'] = {"arg-his crosslink"                        :-5}
mods['crosslinks']['HS'] = {"alaninohistidine crosslink"               :-18}

mods['nn'] =               {"selenocystine"                            :121,
                            "selenomethionine"                         :149,
                            "homoserine"                               :95,
                            "naphthylalanine"                          :197}

mods['S'] = {        "dehydroalanine"                           :-18,
                     "pyruvoyl-serine"                          :-16,
                     "formylglycine"                            :-2,
                     "O-methylation"                            :14,
                     "O acetylation"                            :42,
                     "selenocysteine"                           :64,
                     "dehydroalanine phosphoserine beta elimination":69,
                     "Sulphonation with SO3H"                   :80,        # hydrolysis to ether [wikipedia] is that all?
                     "Phosphorylation"                          :80,
                     "O-GlcNAc-1-phosphorylation"               :266,
                     "O-pantetheinephosphorylation"             :324}

mods['T'] = {        "methyldehydroalanine"                     :-18,
                     "O-methylation"                            :14,
                     "Sulphonation with SO3H"                   :80,
                     "Phosphorylation"                          :80}

mods['D'] = {        "Phosphorylation"                          :80,
                     "Phenyl Ester"                             :76,
                     "Glycerol Ester"                           :74,
                     "transamidation with piperidine"           :67,
                     "beta-methylthiolation ribosomal"          :46,
                     "Potassium"                                :38,
                     "hydroxylation of beta-C"                  :16,
                     "sodium"                                   :22}

mods['Q'] = {        "deamidation"                              :1,
                     "pyroglutamic acid"                        :-17}

mods['N'] = {        "N-methylation"                            :14,
                     "deamidation"                              :1,
                     "succinimide"                              :-17}

mods['W'] = {        "2,4-BisTrp-6,7-dione"                     :28,
                     "reduction of indole double bond"          :2,
                     "oxidation to kynurenine"                  :4,
                     "formaldehyde adduct"                      :12,
                     "hydroxylation of beta-C"                  :16,
                     "6,7 Dione"                                :30,
                     "double oxidation"                         :32}


mods['Y'] = {        "O-(8 alpha-Flavin [FAD])"                 :783,
                     "3,5,3'-triiodothyronine"                  :470,
                     "beta-glycosyloxy"                         :177,
                     "3-Bromination with 79Br"                  :78,
                     "nitro NO2"                                :45,
                     "3,4,6-trihydroxy-Phe"                     :32,
                     "3,4-dihydroxy-phenylalanine"              :16,
                     "3-chlorination with 35Cl"                 :34,
                     "3-chlorination with 37Cl"                 :36,
                     "3,5-dichlorination with 35Cl"             :68,
                     "3,5-dichlorination with 35Cl and 37Cl"    :70,
                     "Sulphation"                               :80,
                     "Phosphorylation"                          :80,
                     "3-Bromination with 81Br"                  :80,
                     "C3 iodination"                            :126,
                     "3,5-dibromination with 79Br"              :156,
                     "3,5-dibromination with 79Br,81Br mixture" :158,
                     "3,5-dibromination with 81Br"              :160,
                     "3,5-diiodination"                         :252,
                     "O-uridinylylation"                        :306}

mods['P'] = {        "hydroxylation of C3 or C4"                :16,
                     "oxidation to gamma-glutamyl semialdehyde" :16,
                     "3,4-dihydroxylation"                      :32,
                     "oxidation to E"                           :32}

mods['H'] = {        "oxohistidine"                             :16,
                     "phosphate/sulphate adduct"                :98,
                     "C4 iodination"                            :126,
                     "N-theta-(ADP-ribosyl) dipthamide"         :648,
                     "N-theta and N-pi-(8-alpha-flavin)"        :784}

mods['F'] = {        "L-o-bromination with 79Br"                :78,
                     "L-o-bromination with 81Br"                :80,
                     "Sulphonation with SO3H"                   :80,
                     "beta-glycosyloxy"                         :177}

mods['G'] = {        "myristoleylation"                         :208,
                     "myristoylation"                           :208}

mods['nterm'] = {    "nterm NHS ester fluorescein label"        :359,
                     "nterm methylation"                        :14,
                     "nterm formylation"                        :28,
                     "nterm acetylation"                        :42,
                     "nterm carbamylation"                      :43,
                     "nterm pyruvate"                           :70,
                     "nterm trifluoroacetyl"                    :96}

mods['cterm'] = {    "cterm amide"                              :-1,
                     "cterm O-methylation"                      :14,
                     "cterm sodium"                             :22,
                     "cterm potassium"                          :38,
                     "cterm phenyl ester"                       :76,
                     "cterm O-ADP-ribosylation"                 :541}

mods['X'] = {        "nonspecific ethylation"                   :28,
                     "nonspecific acetylation"                  :42,
                     "nonspecific piperidine adduct"            :51,
                     "nonspecific t-butyl ester and t-butyl"    :56,
                     "nonspecific sarcosyl detergent"           :71}
import copy
mods['nsp'] = copy.deepcopy(mods['X'])

import pickle
file = open('C:/Users/Dude/Desktop/SPADE/Tools/MassSpectrometry/mods.pkl', 'wb')
pickle.dump(mods, file, 2)
file.close()

