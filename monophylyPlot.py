#!/usr/bin/env python

# Author: Jesper Svedberg (jesper.svedberg@ebc.uu.se)
# Version: 0.1


from ete3 import Tree
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse
import gzip
import itertools
import os

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'svgfont'
mpl.rcParams['ps.fonttype'] = 42

iwh_colors30 = ["#ff0031",
                "#008e08",
                "#7a008d",
                "#edbb04",
                "#0043a9",
                "#f87800",
                "#a76edc",
                "#b1d263",
                "#ef33a8",
                "#5aa85b",
                "#7e2f86",
                "#96933d",
                "#697bc1",
                "#ce8c41",
                "#79a0e7",
                "#c05e39",
                "#01a0c2",
                "#d53f60",
                "#7dd1b3",
                "#923461",
                "#3f723c",
                "#c873b7",
                "#404100",
                "#ffaad5",
                "#875929",
                "#34659e",
                "#fefed0",
                "#832f2b",
                "#5b3f79",
                "#d47a7b"]

iwh_colors15 = ["#ff0031",
                "#008e08",
                "#7a008d",
                "#edbb04",
                "#0043a9",
                "#f87800",
                "#01a0c2",
                "#ef33a8",
                "#b1d263",
                "#697bc1",
                "#c05e39",
                "#7dd1b3",
                "#ffaad5",
                "#404100",
                "#fefed0"]

def avg(lst):
    return sum(lst)/len(lst)

def nlsort(series,categories):
    lst=[]
    ncat = len(categories)
    for n in series:
        found=False
        cat=0
        for i in categories:
            if n <= i:
                lst.append(cat)
                found=True
                break
            cat += 1
        if not found:
            lst.append(np.nan)
            
    return pd.Series(lst)

def swDF(df,data,coord,cats,size=50,step=25):
    cdict={}
    for n in cats:
        cdict[n] = []
        
    ldf=len(df[data])
    wins = [df[data][x:(x+size)].tolist() for x in range(0,ldf,step)]
    for w in wins:
        for c in cats:
            cdict[c].append(w.count(c))
    
    out = pd.DataFrame.from_dict(cdict)
    out.index = [avg(df[coord][x:(x+size)]) for x in range(0,ldf,step)]
    
    return out

def swDFscale(df,data,coord,scale,cats,size=50,step=25):
    cdict={}
    for n in cats:
        cdict[n] = []
        
    ldf=len(df[data])
    wins = [df[data][x:(x+size)].tolist() for x in range(0,ldf,step)]
    swins = [df[scale][x:(x+size)].tolist() for x in range(0,ldf,step)]
    for i,w in enumerate(wins):
        for c in cats:
            cdict[c].append(w.count(c)*(swins[i]/100.0))
    
    out = pd.DataFrame.from_dict(cdict)
    out.index = [avg(df[coord][x:(x+size)]) for x in range(0,ldf,step)]
    
    return out

def comblist(lst):
    out=[]
    combs=reversed(range(2,len(lst)+1))
    for i in combs:
        for c in itertools.combinations(lst, i):
            out.append(list(c))
    return out

def comblist2(lst,start=1):
    out=[]
    combs=reversed(range(start,len(lst)+1))
    for i in combs:
        for c in itertools.combinations(lst, i):
            out.append(list(c))
    return out

def genCombDict(comblist):
    odict = {}
    for cl in comblist:
        nstr = "-".join(cl)
        odict[nstr] = cl
    return odict

def genCombNames(comblist):
    olist = []
    for cl in comblist:
        nstr = "-".join(cl)
        olist.append(nstr)
    return olist


#Args: file handle for tree file, list of combinations to test
def combinMonos(trees, combins):
    names=[]
    for i in combins:
        ostr = "-".join(i)
        names.append(ostr)
    
    tlist=[]
    slist=[]
    for tr in trees:
        ostr=""
        t = Tree(tr)
        cstr=""
        svals=[]
        for i, c in enumerate(combins):
            if t.check_monophyly(values=c, target_attr="name")[0]:

                if len(c) > 2:
                    tlist.append(names[i])
                    ostr="1"
                    slist.append(t.get_common_ancestor(c).support)
                    break
                else:
                    if not cstr:
                        cstr = names[i]
                        svals.append(t.get_common_ancestor(c).support)
                    else:
                        cstr += "-" + names[i]
                        svals.append(t.get_common_ancestor(c).support)

        if ostr=="1":
            continue
        ostr = cstr
        if ostr:
            tlist.append(ostr)
            slist.append(avg(svals))
        else:
            tlist.append("NA")
            slist.append(0)

    return pd.Series(tlist), pd.Series(slist)

def treeify(names,top=1):
    if top==1:
        for i, n in enumerate(names):
            if i == 0:
                ostr = n
            else:
                ostr = "(" + ostr + "," + n + ")"
        ostr += ";"
    elif top==2:
        ostr = "((" + names[0] + "," + names[1] + "),(" + names[2] + "," + names[3] + "));"
    #print ostr
    t = Tree(ostr)
    return t.get_topology_id(), ostr

def getTopos(nlist):
    perms = list(itertools.permutations(nlist))
    odict={}
    if len(nlist) == 3:
        for p in perms:
            a,b = treeify(p)
            odict[a]=b
    elif len(nlist) == 4:
        for p in perms:
            a,b = treeify(p)
            odict[a]=b
        
        firstn=nlist[0]
        for i, n in enumerate(nlist[1:]):
            tl=list(nlist[1:])
            tl.remove(n)            
            p = [firstn, n, tl[0], tl[1]]
            a,b = treeify(p, 2)
            odict[a]=b
    
    return odict

def comblist3(lst,minim=2):
    out=[]
    combs=reversed(range(minim,len(lst)+1))
    for i in combs:
        for c in itertools.combinations(lst, i):
            out.append(list(c))
    return out

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z




def main():
    parser = argparse.ArgumentParser(description='Script for analyzing monophyly of subclade.')
    parser.add_argument('-t','--treefile', help='Gzipped Newick tree file from phyml_sliding_window.py.')
    parser.add_argument('-s','--statsfile', help='Stats file from phyml_sliding_window.py.')
    parser.add_argument('-m','--mono', help='Group to check for monophyly. Names separated by commas, without spaces.')
    parser.add_argument('-e','--extra', help='Extra strains to check if grouping with monophyletic group.')
    parser.add_argument('-o','--outformat', help='Output image format. PNG or SVG. Default: PNG', default='PNG')
    parser.add_argument('-x','--markers', help='List of sites to be marked by vertical lines.')
    parser.add_argument('-w','--window',help='Number of trees in sliding window plots. Default=20', type=int, default=50)
    parser.add_argument('-v','--step',help='Number of trees in sliding window plots. Default=20', type=int, default=25)
    parser.add_argument('-p','--prefix', help='Output prefix. Default: same as statsfile')
    parser.add_argument('-d','--datafile', help='CSV exported from earlier run.')
    parser.add_argument('-c','--varcolor', help='Allow colors for topology plots to change between script runs.',dest='colors', action='store_false', default=False)
    args = parser.parse_args()

    treefile = args.treefile
    statsfile = args.statsfile

    window = args.window
    step = args.step

    
    if args.prefix is None:
        prefix = statsfile
    else:
        prefix = args.prefix
    
    outdir =  prefix + ".output"
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    outformat = args.outformat.lower()


    if args.datafile is None:
        datafile = ""
    else:
        datafile = args.datafile


    monoGroup = args.mono.split(",")
    mCombs = comblist(monoGroup)
    mCombsNames = genCombNames(mCombs)

    if args.extra is None:
        extra = []
    else:
        extra = args.extra.split(",")
    mon_ext = monoGroup + extra
    
    meCombs = []
    for i in mCombs:
        meCombs.append(i + extra)
    meCombsNames = genCombNames(meCombs)

    melCombs=[]
    for ll in comblist2(monoGroup):
        melCombs.append(ll+extra)
    melCombsNames = genCombNames(melCombs)


    mGroups_nc =  comblist3(monoGroup,3)
    mGroups_topos = {}
    for c in mGroups_nc:
        mGroups_topos = merge_two_dicts(mGroups_topos, getTopos(c))


    if args.markers is None:
        markers = []
    else:
        markers = [int(x) for x in args.markers.split(",")]
    
    # Comparable colors
    compColors = args.colors

    if not datafile:
        treeDF=pd.read_csv(statsfile,delimiter="\t")
        treeList = [Tree(tr) for tr in gzip.open(treefile)]
        #treeDF["tree"] = pd.Series(treeList)


        is_mono = []
        is_mon_ext = []
        is_sup = []
        treesup = []
        nodesize = []
        tlist=[]
        slist=[]
        me_tlist=[]
        me_slist=[]
        mel_tlist=[]
        mel_slist=[]
        topo_list=[]
        topo_sup_list=[]

        for xx, tree in enumerate(treeList):
            mono = tree.check_monophyly(values=monoGroup, target_attr="name")[0]
            is_mono.append(mono)
            
            mex = tree.check_monophyly(values=mon_ext, target_attr="name")[0]
            is_mon_ext.append(mex)
            
            sup = []
            for node in tree.traverse():
                if not node.is_leaf() and not node.is_root():
                    sup.append(node.support)

            treesup.append(avg(sup))
            mG_ca = tree.get_common_ancestor(monoGroup)
            is_sup.append(mG_ca.support)
            nodesize.append(len(mG_ca))

            ostr=""
            cstr=""
            svals=[]
            for i, c in enumerate(mCombs):
                if tree.check_monophyly(values=c, target_attr="name")[0]:
                    
                    if len(c) > 2:
                        tlist.append(mCombsNames[i])
                        ostr="1"
                        slist.append(tree.get_common_ancestor(c).support)
                        break
                    else:
                        if not cstr:
                            cstr = mCombsNames[i]
                            svals.append(tree.get_common_ancestor(c).support)
                        else:
                            cstr += "+" + mCombsNames[i]
                            svals.append(tree.get_common_ancestor(c).support)
            
            if ostr != "1":
                ostr = cstr
                if ostr:
                    tlist.append(ostr)
                    slist.append(avg(svals))
                else:
                    tlist.append("NA")
                    slist.append(0)
            

            if tlist[-1] != "NA" and tlist[-1].count("-") > 1 and tlist[-1].count("-") < 4 and tlist[-1].count("+") < 1:
                #print tlist[-1].split("-")
                #print xx
                ca = tree.get_common_ancestor(tlist[-1].split("-"))
                
                sup=[]
                for node in ca.traverse():
                    if not node.is_leaf() and not node.is_root():
                        sup.append(node.support)
                topo_sup_list.append(avg(sup))
                #print ca.get_ascii(attributes=["name","support"], show_internal=True)
                topo_list.append(mGroups_topos[ca.get_topology_id()])

            else:
                topo_sup_list.append(np.nan)
                topo_list.append(np.nan)
            



            
            if extra:
                ostr=""
                cstr=""
                svals=[]
                for i, c in enumerate(meCombs):
                    if tree.check_monophyly(values=c, target_attr="name")[0]:
                        
                        if len(c) > 2:
                            me_tlist.append(meCombsNames[i])
                            ostr="1"
                            me_slist.append(tree.get_common_ancestor(c).support)
                            break
                        else:
                            if not cstr:
                                cstr = meCombsNames[i]
                                svals.append(tree.get_common_ancestor(c).support)
                            else:
                                cstr += "+" + meCombsNames[i]
                                svals.append(tree.get_common_ancestor(c).support)
                
                if ostr != "1":
                    ostr = cstr
                    if ostr:
                        me_tlist.append(ostr)
                        me_slist.append(avg(svals))
                    else:
                        me_tlist.append("NA")
                        me_slist.append(0)
                
                
                ostr=""
                cstr=""
                svals=[]
                for i, c in enumerate(melCombs):
                    if tree.check_monophyly(values=c, target_attr="name")[0]:
                        
                        if len(c) > 2:
                            mel_tlist.append(melCombsNames[i])
                            ostr="1"
                            mel_slist.append(tree.get_common_ancestor(c).support)
                            break
                        else:
                            if not cstr:
                                cstr = melCombsNames[i]
                                svals.append(tree.get_common_ancestor(c).support)
                            else:
                                cstr += "+" + melCombsNames[i]
                                svals.append(tree.get_common_ancestor(c).support)
                
                if ostr != "1":
                    ostr = cstr
                    if ostr:
                        mel_tlist.append(ostr)
                        mel_slist.append(avg(svals))
                    else:
                        mel_tlist.append("NA")
                        mel_slist.append(0)

        
        treeDF["is_monophyletic"] = pd.Series(is_mono)
        treeDF["is_mon_ext"] = pd.Series(is_mon_ext)
        treeDF["mono_support"] = pd.Series(is_sup)
        treeDF["tree_support"] = pd.Series(treesup)
        treeDF["nodesize"] = pd.Series(nodesize)
        treeDF["monocats"] = pd.Series(tlist)
        treeDF["monosupport"] = pd.Series(slist)
        treeDF["monosupport2"] = treeDF.apply(lambda row: row['mono_support'] if (row['monosupport']==0) else row['monosupport'], axis=1)
        treeDF["me_cats"] = pd.Series(me_tlist)
        treeDF["me_support"] = pd.Series(me_slist)
        treeDF["mel_cats"] = pd.Series(mel_tlist)
        treeDF["mel_support"] = pd.Series(mel_slist)
        treeDF["mel_support2"] = treeDF.apply(lambda row: row['mono_support'] if (row['mel_support']==0) else row['mel_support'], axis=1)
        treeDF["sk2rec"] =  pd.Series(topo_list)
        treeDF["sk2rec_sup"] =  pd.Series(topo_sup_list)

        
        outdata = outdir + "/" + prefix + ".out_data.csv"
        treeDF.to_csv(outdata)
    else:
        treeDF = pd.read_csv(datafile, index_col=0)

    startCoord = treeDF.iloc[0].start
    endCoord = treeDF.iloc[-1].end

    #colors = iwh_colors15
    #colors=["#aaaaaa", "#ffc000", "#25f9ef", "#9e06e0", "#ff3030", "#ff00f2", "#002eff", "#5de006", "#fffcad", "#ffd6fa", "#6b0000", "#d0cfd1", "#611870", "#1c8260", "#ff0000", "#7f8400", "#000000", "#84007d", "#ffffff", "#00ff00", "#0000ff", "#ff00ff","#BF3EFF","#4B0082","#B0C4DE","#FFF68F","#8B5A00","#8B0000","#4A708B","#9BCD9B","#FFDAB9","#8B7765","#8B4513","#8E388E"]
    cc=['7402-7426+7427-7429','7402-7427+7426-7429','7402-7429+7426-7427']
    mColDict = dict(zip(mCombsNames+cc+['NA'],iwh_colors15))
    meColDict = dict(zip(meCombsNames+cc+['NA'],iwh_colors15))
    melColDict = dict(zip(melCombsNames+cc+['NA'],iwh_colors15))

    if compColors:
        topos=[]
        for key, val in mGroups_topos.iteritems():
            topos.append(val)
        #topos.sort()
        topoColDict = dict(zip(topos,iwh_colors30))
    else:
        uTopos = treeDF['sk2rec'].unique()
        uTopos = filter(lambda v: v==v, uTopos)
        topoColDict = dict(zip(uTopos,iwh_colors30))

    treeDF[treeDF["is_monophyletic"]==False].plot(kind="scatter", x="mid", y="mono_support", c="nodesize", cmap="plasma", figsize=(20,8))

    if markers:
        for n in markers:
            plt.axvline(n)
    sns.plt.ylim(0, 120)
    sns.plt.xlim(startCoord, endCoord)

    outname = outdir + "/" + prefix + ".01." + outformat
    plt.savefig(outname, bbox_inches='tight', dpi=200)

    #colors=["#aaaaaa", "#ffc000", "#25f9ef", "#c7fcd2", "#ff3030", "#237034", "#002eff", "#5de006", "#fffcad", "#ffd6fa", "#6b0000", "#754f0c", "#611870", "#1c8260"]
    #sns.set_palette(colors)
    #fig, ax2 = plt.subplots(figsize=(15,3))
    SFONT=18*0.783

    rc={'font.size': SFONT, 'axes.labelsize': SFONT, 'legend.fontsize': SFONT, 
        'axes.titlesize': SFONT, 'xtick.labelsize': SFONT, 'ytick.labelsize': SFONT}

    with sns.plotting_context("paper"):
        sns.set(rc=rc)
        #plt.figure(fig_size=(15, 3))
        sns.lmplot(data=treeDF, x="mid", y="monosupport2", hue="monocats", fit_reg=False, legend=False, palette=mColDict, size=4, aspect=4.5)
        for n in markers:
            plt.axvline(n, color='r')
        plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0., fontsize=12)
        sns.plt.ylim(0, 110)
        sns.plt.xlim(startCoord, endCoord)
        plt.ylabel('')
        outname = outdir + "/" + prefix + ".02_SUB." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



        cats = treeDF['monocats'].unique()
        phylDF = swDF(treeDF,'monocats','mid',cats,window,step)
        phylDF.plot(kind='area', stacked=False, figsize=(20,4), color=[mColDict.get(x, '#888888') for x in phylDF.columns])
        #for n in markers:
        #    plt.axvline(n, color='k')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        outname = outdir + "/" + prefix + ".03_SUB." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



        #colors=["#aaaaaa", "#ffc000", "#25f9ef", "#9e06e0", "#ff3030", "#237034", "#002eff", "#5de006", "#fffcad", "#ffd6fa", "#6b0000", "#d0cfd1", "#611870", "#1c8260"]
        #sns.set_palette(colors)
        sns.lmplot(data=treeDF, x="mid", y="me_support", hue="me_cats", fit_reg=False, size=7, aspect=3, legend=False, palette=meColDict)
        for n in markers:
            plt.axvline(n, color='r')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        sns.plt.ylim(0, 110)
        sns.plt.xlim(startCoord, endCoord)
        outname = outdir + "/" + prefix + ".04." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



        cats = treeDF['me_cats'].unique()
        #colors=["#aae2d2", "#ffc000", "#25f9ef", "#9e06e0", "#ff3030", "#237034", "#002eff", "#5de006", "#fffcad", "#d0cfd1", "#6b0000", "#d0cfd1", "#611870", "#1c8260"]
        #sns.set_palette(colors)
        phylDFme = swDF(treeDF,'me_cats','mid',cats,window,step)
        phylDFme.plot(kind='area', stacked=False, figsize=(16,5), color=[meColDict.get(x, '#888888') for x in phylDFme.columns])
        for n in markers:
            plt.axvline(n, color='k')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        outname = outdir + "/" + prefix + ".05." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



        #colors=["#aaaaaa", "#ffc000", "#25f9ef", "#9e06e0", "#ff3030", "#237034", "#002eff", "#5de006", "#fffcad", "#ffd6fa", "#6b0000", "#d0cfd1", "#611870", "#1c8260", "#ff0000", "#ffffff", "#000000"]
        #sns.set_palette(colors)
        sns.lmplot(data=treeDF, x="mid", y="mel_support2", hue="mel_cats", fit_reg=False, size=4, aspect=4.5, legend=False, palette=melColDict)
        for n in markers:
            plt.axvline(n, color='r')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0., fontsize=12)
        sns.plt.ylim(0, 110)
        sns.plt.xlim(startCoord, endCoord)
        outname = outdir + "/" + prefix + ".06." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



        #cp=["#aaaaaa", "#ffc000", "#25f9ef", "#9e06e0", "#ff3030", "#ff00f2", "#002eff", "#5de006", "#fffcad", "#ffd6fa", "#6b0000", "#d0cfd1", "#611870", "#1c8260", "#ff0000", "#7f8400", "#000000", "#84007d", "#ffffff"]
        #sns.set_palette(cp)
        sns.lmplot(data=treeDF, x="mid", y="sk2rec_sup", hue="sk2rec", fit_reg=False, size=7, aspect=3, legend=False, palette=topoColDict)
        for n in markers:
            plt.axvline(n, color='r')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        sns.plt.ylim(0, 110)
        sns.plt.xlim(startCoord, endCoord)
        outname = outdir + "/" + prefix + ".07." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



        cats = treeDF['sk2rec'].unique()
        phylDFr = swDF(treeDF,'sk2rec','mid',cats,window,step)
        phylDFr.plot(kind='area', stacked=False, figsize=(16,5) , color=[topoColDict.get(x, '#888888') for x in phylDFr.columns])
        for n in markers:
            plt.axvline(n, color='k')
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        outname = outdir + "/" + prefix + ".08." + outformat
        plt.savefig(outname, bbox_inches='tight', dpi=200)



    # SUBPLOTS of 02 and 03
'''
    fig, (ax1, ax2) = plt.subplots(2, figsize=(22,12))

    sns.lmplot(data=treeDF, x="mid", y="monosupport2", hue="monocats", fit_reg=False, size=7, aspect=2, legend=False, palette=mColDict)
    for n in markers:
        plt.axvline(n, color='r')
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0.)
    
    phylDF.plot(kind='area', stacked=False, figsize=(16,8), color=[mColDict.get(x, '#888888') for x in phylDF.columns])
    for n in markers:
        plt.axvline(n, color='k')
    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    outname = outdir + "/" + prefix + ".subplots." + outformat
    plt.savefig(outname, bbox_inches='tight', dpi=200)'''






    


if __name__ == '__main__':
    main()
