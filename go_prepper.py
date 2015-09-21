def prepper():
    """
    Function that performs the go.obo parser
    Was adjusted to use the version format 1.2 dated 2015-03-27
    """
    go_dict = {}
    with open("go.obo") as f:
        content = f.readlines()
        s = False
        d = {}
        a = 1
    for i in content:  
        
        if (s) and (i != "[Term]\n") and (i != "\n"):
            k = (i.split(": ")[0])
            c = (i.split(": ")[1])
            if d.has_key(k): 
                d[k+str(a)] = c
                a+=1
            else:
                d[k] = c
                
        if (i == "[Term]\n"):
            d = {}
            s = True

        elif (i == "\n"):
            s = False
            a = 1
            if len(d) != 0 :
                gt = d["id"].split("\n")[0]
                go_dict[gt] = d
    return go_dict


def slimmer (go_dict,gset):
    """
    making a goslim dictionary
    """
    go_slims = {}
    for i in go_dict.keys():
        for j in go_dict[i].keys():
            if (j[:6]=="subset") and (go_dict[i][j] == 'goslim_generic\n'):
                if gset.has_key(i):
                    go_slims[i] = gset[i]
                else:
                    print go_dict[i]["name"]

    return go_slims
    


def geneset_association():
    """
    function to verify if the gene is already in the set
    """
    gset = {}
    f = open ("gene_association.goa_human")
    for line in f:
        if line[0] != "!": 
            i = (line.split ('\t'))
            go = i[4]
            gene = []    
            gene.append(i[2])
            gene = set(gene)
        
            if not gset.has_key (go): 
                gset[go] = set()

            if (gene not in gset[go]):
                gset[go] |= gene

    return gset


def traceback (parent,genes,gset,go_dict):
    """
    Function that does a traceback and adds the is_a hierarchy to gene sets
    """
    if not gset.has_key(parent):
        gset[parent] = set()

    gset[parent] |= genes
    pnames = go_dict[parent]
    
    for i in pnames.keys():
        if i[:4] == "is_a":
            child = pnames[i].split(" ! ")[0]
            traceback(child,genes,gset,go_dict)

    
def prepGOs():
    """
    making the GO term dictionnary go_dict
    go_dict will contain keys as GO terms.
    to each key it's dictionnary ("is_a" are stored here!)
    """

    gset = geneset_association()            
    go_dict = prepper()
    go_slims = slimmer(go_dict,gset)

    for j in gset.keys():
        traceback (j,gset[j],gset,go_dict)

    pickle.dump(go_slims,open("goslim.p","wb"))
    pickle.dump(gset,open("Gset.p","wb"))
    pickle.dump(go_dict, open ("go_dict.p","wb"))

