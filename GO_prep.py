import cPickle as pickle
import os.path
from GOprep import prepper,slimmer,geneset_association,traceback



#making the GO term dictionnary go_dict
#go_dict will contain keys as GO terms.
#to each key it's dictionnary ("is_a" are stored here!)


if os.path.exists("go_dict.p") and os.path.exits("Gset.p") and os.path.exists("goslim.p"):
    #checking if the pickle file exists already
    go_dict = pickle.load(open("go_dict.p"))


else:
    gset = geneset_association()            
    go_dict = prepper()
    go_slims = slimmer(go_dict,gset)

    for j in gset.keys():
        traceback (j,gset[j],gset,go_dict)

    pickle.dump(go_slims,open("goslim.p","wb"))
    pickle.dump(gset,open("Gset.p","wb"))
    pickle.dump(go_dict, open ("go_dict.p","wb"))

    
            

# # tests, nothing to see here, move along
# test = {}
# test["GO:0019992"] = gset["GO:0019992"][:1]
# test["GO:0008289"] = gset["GO:0008289"][:1]
# test["GO:0003674"] = gset["GO:0003674"][:1]


            

    
