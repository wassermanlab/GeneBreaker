# n.b. the base_pb2 and phenopackets_pb2 files are located in the target/generated-sources/protobuf/python
# directory once the project has built. This file *should* run when done so via IntelliJ/PyCharm but otherwise
# some file hackery will be required.

import os,sys

from google.protobuf.json_format import Parse, MessageToJson
from google.protobuf.timestamp_pb2 import Timestamp

from base_pb2 import *
from phenopackets_pb2 import *

# Main - we're going to create a simple phenopacket, write it out as JSON then read it back in as an object
subject = Individual(id="Zaphod", sex="MALE", date_of_birth=Timestamp(seconds=-123456798))
phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HG2G:00001", label="Hoopy")),
                       PhenotypicFeature(type=OntologyClass(id="HG2G:00002", label="Frood"))
                       ]
phenopacket = Phenopacket(id="PPKT:1", subject=subject, phenotypic_features=phenotypic_features)

test_json_file = "test.json"
with open(test_json_file, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)
#
#with open(test_json_file, 'r') as jsfile:
#    round_trip = Parse(message=Phenopacket(), text=jsfile.read())
#    print(round_trip)
#
#    individual = round_trip.subject
#    print("Meet", individual.id)
#
#    if individual.sex:
#        sex = Sex.Name(individual.sex)
#        if sex == 'MALE':
#            print("He is", sex)
#        elif sex == 'FEMALE':
#            print("She is", sex)
#
#    print("{}'s friends describe him as a:".format(individual.id))
#    for phenotypic_feature in round_trip.phenotypic_features:
#        term = phenotypic_feature.type
#        print("{} [{}]".format(term.label, term.id))
#
# For GeneBreaker Scenarios
# Case 1
# Male, 2 years old, WAS, clinvar:424355, inheritance: x-linked hemizygous recessive, 
# hpo terms: 'HP:0000964', 'HP:0001873', 'HP:0000388', 'HP:0002248'
# hpo mappings: "Eczema,Thrombocytopenia,Otitis media, Hematemesis"
subject = Individual(id="Patient1", sex="MALE")
# age not recognized???

# doesn't work for some reason...
#PhenotypicFeature(type=description("A two year old boy presents to his pediatrician with bloody diarrhea, eczema, petechiae and easy bruising with no associated trauma. He has had multiple treatments for otitis media (ear infections). Investigations reveal thrombocytopenia (decreased number of platelets) and abnormal morphology (small platelets). The proband is an only child. Both parents do not have any significant medical problems. The family history is significant for a maternal uncle with autoimmune problems and a recent diagnosis of lymphoma."))]

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0000964", label="Eczema")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001873", label="Thrombocytopenia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000388", label="Otitis media")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002248", label="Hematemesis"))]
phenopacket = Phenopacket(id="Case1",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case1_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 2
# Male, 41 year old, MSH2, clinvar:91055, inheritance: autosomal dmoninant de novo
# hpo terms: HP:0200008', 'HP:0001250', 'HP:0001276', 'HP:0012378', 'HP:0002027'
# hpo mappings: Intestinal polyposis, Seizure, Hypertonia, Fatigue, Abdominal pain
subject = Individual(id="Patient2", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0200008", label="Intestinal polyposis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001250", label="Seizure")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001276", label="Hypertonia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0012378", label="Fatigue")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002027", label="Abdominal pain"))]
phenopacket = Phenopacket(id="Case2",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case2_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)

# Case 3
# Female, 14 years old, MALT1, clinvar:662739, inheritance: autosomal recessive homozygous
# hpo terms HP:0000964', 'HP:0001047', 'HP:0001581', 'HP:0004386', 'HP:0002090', 'HP:0002205'
# hpo mappings: Eczema, Atopic dermatitis, Recurrent skin infections, Gastrointestinal inflammation, Pneumonia, Recurrent respiratory infections
subject = Individual(id="Patient3", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0000964", label="Eczema")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001047", label="Atopic dermatitis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001581", label="Recurrent skin infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0004386", label="Gastrointestinal inflammation")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002205", label="Recurrent respiratory infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002090", label="Pneumonia"))]
phenopacket = Phenopacket(id="Case3",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case3_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 4
# Male, 2 years old, CFTR, clinvar: 618298;554293, inheritance: autosomal recessive compound heterozygous
# hpo terms: HP:0002613', 'HP:0002721', 'HP:0002024', 'HP:0002206', 'HP:0002205', 'HP:0001738'
# hpo mappings: Biliary cirrhosis, Immunodeficiency, Malabsorption, Pulmonary fibrosis, Recurrent respiratory infections, Exocrine pancreatic insufficiency
subject = Individual(id="Patient4", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0002613", label="Biliary cirrhosis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002721", label="Immunodeficiency")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002024", label="Malabsorption")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002206", label="Pulmonary fibrosis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002205", label="Recurrent respiratory infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001738", label="Exocrine pancreatic insufficiency"))]
phenopacket = Phenopacket(id="Case4",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case4_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 5
# Female, 2 years old, MECP2, clinvar: 424171, inheritance: x-linked dominant de novo
# hpo terms: HP:0001250', 'HP:0001257', 'HP:0005484', 'HP:0002187'
# hpo mappings: Seizure, Spasticity, Postnatal microcephaly, Intellectual disability profound
subject = Individual(id="Patient5", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001250", label="Seizure")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001257", label="Spasticity")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0005484", label="Postnatal microcephaly")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002187", label="Intellectual disability profound"))]
phenopacket = Phenopacket(id="Case5",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case5_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 6
# Male, 6 years old, ABCD1, clinvar: 458629, inheritance: x-linked recessive hemizygous de novo
# hpo terms: 'HP:0001250', 'HP:0000709', 'HP:0008207', 'HP:0002180', 'HP:0002500'
# hpo mappings: Mental deterioration, Seizure, Psychosis, Primary adrenal insufficiency, Neurodegeneration, Abnormality of the cerebral white matter
subject = Individual(id="Patient6", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001250", label="Seizure")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002180", label="Psychosis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000709", label="Primary adrenal insufficiency")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0008207", label="Neurodegeneration")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002500", label="Abnormality of the cerebral white matter"))]
phenopacket = Phenopacket(id="Case6",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case6_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 7
# Female, EP300, clinvar 666310, inheritance: autosomal dominant de novo
# ‘HP:002342’, ‘HP:0002553’, ‘HP:0000494’, ‘HP:0000448’, ‘HP:0000347’, ‘HP:0011304’
# Moderate intellectual disability, arched eyebrows, downslanting palpebral fissures, prominent nose, micrognathia, broad thumbs. Father with significant learning problems and similar facial features.
subject = Individual(id="Patient7", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:002342", label="Moderate intellectual disability")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002553", label="arched eyebrows")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000494", label="downslanting palpebral fissures")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000448", label="prominent nose")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000347", label="micrognathia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0011304", label="broad thumbs"))]
phenopacket = Phenopacket(id="Case7",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case7_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 8
# Male, PTPN11, Noonan Syndrome, 40488
# HP:0004322’, ‘HP:0032318’, 'HP:0000475', 'HP:0000368', 'HP:0000316', 'HP:0000445', 'HP:0003196'
# Short stature, congenital heart disease, wide neck, low set posteriorly rotated ears, widely spaced eyes, short and broad nose. Negative family history.
subject = Individual(id="Patient8", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0004322", label="Short stature")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0032318", label="congenital heart disease")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000475", label="wide neck")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000368", label="low set posteriorly rotated ears")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000316", label="widely spaced eyes")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000445", label="broad nose")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0003196", label="short nose"))]
phenopacket = Phenopacket(id="Case8",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case8_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Case 9
# Female, COL2A1, 449397
# ‘HP:0011800’, ‘HP:0008625’, ‘HP:0000545’, ‘HP:0000175’,
# Flat midface, sensorineural hearing loss, myopia, cleft palate. Father with severe myopia. 
subject = Individual(id="Patient9", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0011800", label="Flat midface")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0008625", label="sensorineural hearing loss")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000545", label="myopia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000175", label="cleft palate"))]
phenopacket = Phenopacket(id="Case9",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case9_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)



# Case 10
# Male, ANKRD11, 658724
# ‘HP:0001249’, ‘HP:0004322’, ‘HP:0000252’, ‘HP:0000325’, ‘HP:0000400’, ‘HP:0000316’, ‘HP:0000637’, ‘HP:0001572’
# Intellectual disability, short stature, microcephaly, triangular face, large and prominent ears, hypertelorism, long palpebral fissures, macrodontia. Mother similarly affected.
subject = Individual(id="Patient10", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001249", label="Intellectual disability")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0004322", label="short stature")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000252", label="microcephaly")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000325", label="triangular face")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000400", label="large and prominent ears")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000316", label="hypertelorism")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000637", label="long palpebral fissures")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001572", label="macrodontia"))]
phenopacket = Phenopacket(id="Case10",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case10_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# Now for the inheritance testing cases
# JAK1
# HP:0000964; HP:0001047; HP:0001880; HP:0001508; HP:0032064
subject = Individual(id="JAK1", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0000964", label="Eczema")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001047", label="Atopic dermatitis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001880", label="Eosinophilia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001508", label="Failure to thrive")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0032064", label="Gastrointestinal eosinophilia"))]
phenopacket = Phenopacket(id="JAK1",subject=subject, phenotypic_features=phenotypic_features)

json_out = "JAK1_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# MSH2
# hpo terms: HP:0200008', 'HP:0001250', 'HP:0001276', 'HP:0012378', 'HP:0002027'
# hpo mappings: Intestinal polyposis, Seizure, Hypertonia, Fatigue, Abdominal pain
subject = Individual(id="MSH2", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0200008", label="Intestinal polyposis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001250", label="Seizure")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001276", label="Hypertonia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0012378", label="Fatigue")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002027", label="Abdominal pain"))]
phenopacket = Phenopacket(id="MSH2",subject=subject, phenotypic_features=phenotypic_features)

json_out = "MSH2_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)





# MALT1
# hpo terms HP:0000964', 'HP:0001047', 'HP:0001581', 'HP:0004386', 'HP:0002090', 'HP:0002205'
# hpo mappings: Eczema, Atopic dermatitis, Recurrent skin infections, Gastrointestinal inflammation, Pneumonia, Recurrent respiratory infections
subject = Individual(id="MALT1", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0000964", label="Eczema")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001047", label="Atopic dermatitis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001581", label="Recurrent skin infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0004386", label="Gastrointestinal inflammation")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002205", label="Recurrent respiratory infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002090", label="Pneumonia"))]
phenopacket = Phenopacket(id="MALT1",subject=subject, phenotypic_features=phenotypic_features)

json_out = "MALT1_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)




# CFTR
# hpo terms HP:0000964', 'HP:0001047', 'HP:0001581', 'HP:0004386', 'HP:0002090', 'HP:0002205'
# hpo mappings: Eczema, Atopic dermatitis, Recurrent skin infections, Gastrointestinal inflammation, Pneumonia, Recurrent respiratory infections
subject = Individual(id="CFTR", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0000964", label="Eczema")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001047", label="Atopic dermatitis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001581", label="Recurrent skin infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0004386", label="Gastrointestinal inflammation")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002205", label="Recurrent respiratory infections")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002090", label="Pneumonia"))]
phenopacket = Phenopacket(id="CFTR",subject=subject, phenotypic_features=phenotypic_features)

json_out = "Case3_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)



# INPP5E
# HP:0001252', 'HP:0001251', 'HP:0001263', 'HP:0002553', 'HP:0001320', 'HP:0002793'

subject = Individual(id="INPP5E", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001252", label="Muscular hypotonia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001251", label="Ataxia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001263", label="Global developmental delay")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002553", label="Highly arched eyebrow")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001320", label="Cerebellar vermis hypoplasia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002793", label="Abnormal pattern of respiration"))]
phenopacket = Phenopacket(id="INPP5E",subject=subject, phenotypic_features=phenotypic_features)

json_out = "INPP5E_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# MECP2
# hpo terms: HP:0001250', 'HP:0001257', 'HP:0005484', 'HP:0002187'
# hpo mappings: Seizure, Spasticity, Postnatal microcephaly, Intellectual disability profound
subject = Individual(id="MECP2", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001250", label="Seizure")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001257", label="Spasticity")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0005484", label="Postnatal microcephaly")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002187", label="Intellectual disability profound"))]
phenopacket = Phenopacket(id="MECP2",subject=subject, phenotypic_features=phenotypic_features)

json_out = "MECP2_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# WAS
# hpo terms: 'HP:0000964', 'HP:0001873', 'HP:0000388', 'HP:0002248'
# hpo mappings: "Eczema,Thrombocytopenia,Otitis media, Hematemesis"
subject = Individual(id="WAS", sex="MALE")
phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0000964", label="Eczema")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001873", label="Thrombocytopenia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000388", label="Otitis media")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002248", label="Hematemesis"))]
phenopacket = Phenopacket(id="WAS",subject=subject, phenotypic_features=phenotypic_features)

json_out = "WAS_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)



# SLC6A8
# HP:0001290', 'HP:0001270', 'HP:0000252', 'HP:0000718', 'HP:0008583', 'HP:0000540', 'HP:0008583'

subject = Individual(id="SLC6A8", sex="FEMALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001290", label="Generalized hypotonia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0001270", label="Motor delay")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000252", label="Microcephaly")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000718", label="Aggressive behavior")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0008583", label="Underfolded superior helices")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000540", label="Hypermetropia"))]
phenopacket = Phenopacket(id="SLC6A8",subject=subject, phenotypic_features=phenotypic_features)

json_out = "SLC6A8_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)

# ABCD1
# hpo terms: 'HP:0001250', 'HP:0000709', 'HP:0008207', 'HP:0002180', 'HP:0002500'
# hpo mappings: Mental deterioration, Seizure, Psychosis, Primary adrenal insufficiency, Neurodegeneration, Abnormality of the cerebral white matter
subject = Individual(id="ABCD1", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0001250", label="Seizure")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002180", label="Psychosis")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000709", label="Primary adrenal insufficiency")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0008207", label="Neurodegeneration")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0002500", label="Abnormality of the cerebral white matter"))]
phenopacket = Phenopacket(id="Case6",subject=subject, phenotypic_features=phenotypic_features)

json_out = "ABCD1_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)


# SRY
# HP:0012245', 'HP:0011969', 'HP:0000032', 'HP:0000098'
subject = Individual(id="SRY", sex="MALE")

phenotypic_features = [PhenotypicFeature(type=OntologyClass(id="HP:0012245", label="Sex reversal")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0011969", label="Elevated circualting luteinizing hormone level")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000032", label="Abnormality of male external genitalia")),
                        PhenotypicFeature(type=OntologyClass(id="HP:0000098", label="Tall stature"))]
phenopacket = Phenopacket(id="SRY",subject=subject, phenotypic_features=phenotypic_features)

json_out = "SRY_Phenopacket.json"
with open(json_out, 'w') as jsfile:
    json = MessageToJson(phenopacket)
    jsfile.write(json)






sys.exit()

try:
    os.remove(test_json_file)
except OSError as e:
    print("Error: %s - %s." % (e.filename, e.strerror))
