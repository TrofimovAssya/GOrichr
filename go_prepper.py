###go.obo parser

def prepper():
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

#making a goslim dictionary
def slimmer (go_dict,gset):
    go_slims = {}
    for i in go_dict.keys():
        for j in go_dict[i].keys():
            if (j[:6]=="subset") and (go_dict[i][j] == 'goslim_generic\n'):
                if gset.has_key(i):
                    go_slims[i] = gset[i]
                else:
                    print go_dict[i]["name"]

    return go_slims
    
#putting together genesets


#function to verify if the gene is already in the set
def geneset_association():
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

### Use with GO_prep.py
### Function that does a traceback and adds the is_a hierarchy to gene sets

def traceback (parent,genes,gset,go_dict):

    if not gset.has_key(parent):
        gset[parent] = set()

        
    gset[parent] |= genes
    
        
    pnames = go_dict[parent]
    
    for i in pnames.keys():
        if i[:4] == "is_a":
            child = pnames[i].split(" ! ")[0]
            traceback(child,genes,gset,go_dict)

    
