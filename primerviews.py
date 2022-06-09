#coding:utf-8
from django.template import RequestContext,loader
from django.shortcuts import render_to_response
from django.http import HttpResponse,HttpResponseRedirect
from django.views.decorators.csrf import csrf_protect,csrf_exempt
import time
import primer3
import Cookie
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd

localip='www.kunkundashen.cn'
thisyear=time.strftime('%Y',time.localtime(time.time()))
thismonth=time.strftime('%m',time.localtime(time.time()))

def primerdesigner(request):
    return render_to_response('primer3/primer.html',{
                                    'local':localip,
                                    'thisyear':thisyear
        },context_instance=RequestContext(request))
##############导入设计序列################
@csrf_protect  
def seqview(request):
     if request.method=="POST":
        Entrez.email="lxk@kodiabio.com"
        SeqId=request.POST['seqid']
        Start=int(request.POST['start'])
        End=int(request.POST['end'])
        handle=Entrez.efetch(db="nucleotide",rettype="gb",retmote="text",id=SeqId)
        record=SeqIO.read(handle,"gb")
        time.sleep(10)
        handle.close()
        cookie=Cookie.SimpleCookie()
        response=render_to_response('primer3/parser.html',{
                                                       'local':localip,
                                                       'thisyear':thisyear,
                                                       'filetype':'genbank',
                                                       'accessid':record.id,
                                                       'sequence':str(record.seq[Start-1:End]),
                                                       'description':record.description,
                                                       'name':record.name,
                                                       #'dbxrefs':nr.dbxrefs[0],
                                                       'source':record.annotations['source'],
                                                       'organism':record.annotations['organism'],
                                                       'taxonomy':record.annotations['taxonomy'],
                                                       'topology':record.annotations['topology'],
                                                           },context_instance=RequestContext(request))
#         response.set_cookie("seq",replace_RAGTC(str(record.seq[Start-1:])))
        response.set_cookie("des",record.description)
        response.set_cookie("start-end",str(Start)+"-"+str(End))
        response.set_cookie("seqid",SeqId)
        return response
#########远程访问Entrez数据库#####################    
def checkAccessionLen(request,accessionid):
    Entrez.email='lxk@kodiabio.com'
    handle=Entrez.efetch(db='nucleotide',rettype='gb',retmote='text',id=accessionid)
    record=SeqIO.read(handle, 'gb')
    Length=len(record.seq)
    handle.close()
    return HttpResponse(Length)
###################################################################################
####################修改碱基信息###################################################
@csrf_protect
def basemodify(request):
    if request.method=="POST":
        global local
        return render_to_response('primer3/basemodify.html',{
                                                      'local':localip,
                                                      'thisyear':thisyear,
                                                      'sequence':list(request.POST['seq']),
                                                      'description':request.COOKIES.get('des',''),
                                                      },context_instance=RequestContext(request))
        
####################序列信息输入###################################################
@csrf_protect
def seqinput(req):
     if req.method=='POST':
         seqtxt=''.join(req.POST.getlist('seqtxt'))
         seq_args = {
                        'SEQUENCE_ID': 'TARGET_SEGMENT',
                        'SEQUENCE_TEMPLATE': str(seqtxt),
                        'SEQUENCE_INCLUDED_REGION': [0,len(seqtxt)],
                        }
         global_args = {
                        'PRIMER_OPT_SIZE': int(req.POST['PRIMER_OPT_SIZE']),
                        'PRIMER_PICK_INTERNAL_OLIGO': int(req.POST['PRIMER_PICK_INTERNAL_OLIGO']),
                        'PRIMER_INTERNAL_MAX_SELF_END': int(req.POST['PRIMER_INTERNAL_MAX_SELF_END']),
                        'PRIMER_MIN_SIZE': int(req.POST['PRIMER_MIN_SIZE']),
                        'PRIMER_MAX_SIZE': int(req.POST['PRIMER_MAX_SIZE']),
                        'PRIMER_OPT_TM': float(req.POST['PRIMER_OPT_TM']),
                        'PRIMER_MIN_TM': float(req.POST['PRIMER_MIN_TM']),
                        'PRIMER_MAX_TM': float(req.POST['PRIMER_MAX_TM']),
                        'PRIMER_MIN_GC': float(req.POST['PRIMER_MIN_GC']),
                        'PRIMER_MAX_GC': float(req.POST['PRIMER_MAX_GC']),
                        'PRIMER_MAX_POLY_X': int(req.POST['PRIMER_MAX_POLY_X']),
                        'PRIMER_INTERNAL_MAX_POLY_X': int(req.POST['PRIMER_INTERNAL_MAX_POLY_X']),
                        'PRIMER_SALT_MONOVALENT': float(req.POST['PRIMER_SALT_MONOVALENT']),
                        'PRIMER_DNA_CONC': float(req.POST['PRIMER_DNA_CONC']),
                        'PRIMER_MAX_NS_ACCEPTED': int(req.POST['PRIMER_MAX_NS_ACCEPTED']),
                        'PRIMER_MAX_SELF_ANY': int(req.POST['PRIMER_MAX_SELF_ANY']),
                        'PRIMER_MAX_SELF_END': int(req.POST['PRIMER_MAX_SELF_END']),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(req.POST['PRIMER_PAIR_MAX_COMPL_ANY']),
                        'PRIMER_PAIR_MAX_COMPL_END': int(req.POST['PRIMER_PAIR_MAX_COMPL_END']),
                        'PRIMER_PRODUCT_SIZE_RANGE': [[75,200]],
                        }
         primer3_primary_result = primer3.bindings.designPrimers(seq_args, global_args)
         #第一步：将结果转换成以“引物信息”为键值的新形式的“字典”：
         primer3_result_table_dict = {}
         for i in range(primer3_primary_result["PRIMER_PAIR_NUM_RETURNED"]):
             primer_id = str(i)
             for key in primer3_primary_result:
                 if primer_id in key:
                     info_tag = key.replace("_" + primer_id, "")# 要将每个信息中的数字和下划线去掉
                     try:# 就是把不同的引物对结果归到一起
                         primer3_result_table_dict[info_tag]
                     except:
                         primer3_result_table_dict[info_tag] = []
                     finally:
                         if isinstance(primer3_primary_result[key],float):
                             primer3_primary_result[key]=float('%.2f' % primer3_primary_result[key])
                             primer3_result_table_dict[info_tag].append(primer3_primary_result[key])
                         else:
                             primer3_result_table_dict[info_tag].append(primer3_primary_result[key])
         if req.POST['PRIMER_PICK_INTERNAL_OLIGO'] == "0":
             dict_intnl_seqs={"PRIMER_INTERNAL_SEQUENCE":['N/A','N/A','N/A','N/A','N/A']}
         #挑出引物序列    
         for k,v in primer3_result_table_dict.items():
            # dict_intnl_seqs.setdefault("PRIMER_INTERNAL_SEQUENCE", default=['','','','','']) 
             if k=="PRIMER_LEFT_SEQUENCE":
                 dict_left_seqs={"PRIMER_LEFT_SEQUENCE":v}
             if k=="PRIMER_RIGHT_SEQUENCE":
                 dict_right_seqs={"PRIMER_RIGHT_SEQUENCE":v}
                 vv_r_c=[]
                 for vv in v:
                     vv_r_c.append(Seq(vv,IUPAC.unambiguous_dna).reverse_complement())
                 dict_right_RC_seqs={"PRIMER_RIGHT_RC_SEQUENCE":vv_r_c}
             if k=="PRIMER_INTERNAL_SEQUENCE":
                 dict_intnl_seqs={"PRIMER_INTERNAL_SEQUENCE":v}
             if k=="PRIMER_PAIR_PRODUCT_SIZE":
                 dict_prdct_sizes={"PRIMER_PAIR_PRODUCT_SIZE":v}
             if k=="PRIMER_PAIR_COMPL_ANY_TH":
                 dict_compl_anyths={'PRIMER_PAIR_COMPL_ANY_TH':v}
             if k=="PRIMER_PAIR_PENALTY":
                 dict_pair_penaltys={'PRIMER_PAIR_PENALTY':v}

         return render_to_response('primer3/showprimers.html',{
                                                     'local':localip,
                                                     'thisyear':thisyear,
                                                     'seqtxt':str(seqtxt),
                                                     'primer3_primary_result':primer3_result_table_dict,
                                                     'dict_left_seqs':dict_left_seqs,
                                                     'dict_right_seqs':dict_right_seqs,
                                                     'dict_right_RC_seqs':dict_right_RC_seqs,
                                                     'dict_intnl_seqs':dict_intnl_seqs,
                                                     'dict_prdct_sizes':dict_prdct_sizes,
                                                     'dict_compl_anyths':dict_compl_anyths,
                                                     'dict_pair_penaltys':dict_pair_penaltys,
                                                     'description':req.COOKIES.get('des',''),
                                                     'se':req.COOKIES.get('start-end',''),
                                                     'seqid':req.COOKIES.get('seqid',''),
                                                     },context_instance=RequestContext(req))  
##############################################################################################################################
####################序列读取###################################################
@csrf_protect
def seqread(req):
    if req.method=="POST":
        gene_sybol=req.POST['gene_sybol']
        gene_sequence=req.POST['gene_sequence']
        response=render_to_response('primer3/reader.html',{
                                                     'local':localip,
                                                     'thisyear':thisyear,
                                                     'gene_sybol':gene_sybol,
                                                     'gene_sequence':gene_sequence
                                                     },context_instance=RequestContext(req))
        response.set_cookie("des",gene_sybol)
        return response
    
def pcrprincp(req):
    return render_to_response('primer3/pcrprincple.html',{
                                    'local':localip,
                                    'thisyear':thisyear
        },context_instance=RequestContext(req))

def manual(req):
    return render_to_response('primer3/manual.html',{
                                    'local':localip,
                                    'thisyear':thisyear
        },context_instance=RequestContext(req))
def multiprimer(req):
    return render_to_response('multiprimer/multiprimer.html',{
                                    'local':localip,
                                    'thisyear':thisyear
        },context_instance=RequestContext(req))
@csrf_protect
def multiseqview(req):
    if req.method=="POST":
        inputseqs=req.POST.getlist('seqtxt',[])
        count=0
        list_dict_left_seqs=[]
        list_dict_right_seqs=[]
        list_dict_intnl_seqs=[]
        list_dict_prdct_sizes=[]
        list_dict_compl_anyths=[]
        list_dict_pair_penaltys=[]
        list_dict_right_seqs_rc=[] 
        for inputseq in inputseqs:
            seq_args = {
                        'SEQUENCE_ID': 'TARGET_SEGMENT',
                        'SEQUENCE_TEMPLATE': str(inputseq),
                        'SEQUENCE_INCLUDED_REGION': [0,len(inputseq)],
                        }
            global_args = {
                        'PRIMER_OPT_SIZE': int(req.POST['PRIMER_OPT_SIZE']),
                        'PRIMER_PICK_INTERNAL_OLIGO': int(req.POST['PRIMER_PICK_INTERNAL_OLIGO']),
                        'PRIMER_INTERNAL_MAX_SELF_END': int(req.POST['PRIMER_INTERNAL_MAX_SELF_END']),
                        'PRIMER_MIN_SIZE': int(req.POST['PRIMER_MIN_SIZE']),
                        'PRIMER_MAX_SIZE': int(req.POST['PRIMER_MAX_SIZE']),
                        'PRIMER_OPT_TM': float(req.POST['PRIMER_OPT_TM']),
                        'PRIMER_MIN_TM': float(req.POST['PRIMER_MIN_TM']),
                        'PRIMER_MAX_TM': float(req.POST['PRIMER_MAX_TM']),
                        'PRIMER_MIN_GC': float(req.POST['PRIMER_MIN_GC']),
                        'PRIMER_MAX_GC': float(req.POST['PRIMER_MAX_GC']),
                        'PRIMER_MAX_POLY_X': int(req.POST['PRIMER_MAX_POLY_X']),
                        'PRIMER_INTERNAL_MAX_POLY_X': int(req.POST['PRIMER_INTERNAL_MAX_POLY_X']),
                        'PRIMER_SALT_MONOVALENT': float(req.POST['PRIMER_SALT_MONOVALENT']),
                        'PRIMER_DNA_CONC': float(req.POST['PRIMER_DNA_CONC']),
                        'PRIMER_MAX_NS_ACCEPTED': int(req.POST['PRIMER_MAX_NS_ACCEPTED']),
                        'PRIMER_MAX_SELF_ANY': int(req.POST['PRIMER_MAX_SELF_ANY']),
                        'PRIMER_MAX_SELF_END': int(req.POST['PRIMER_MAX_SELF_END']),
                        'PRIMER_PAIR_MAX_COMPL_ANY': int(req.POST['PRIMER_PAIR_MAX_COMPL_ANY']),
                        'PRIMER_PAIR_MAX_COMPL_END': int(req.POST['PRIMER_PAIR_MAX_COMPL_END']),
                        'PRIMER_PRODUCT_SIZE_RANGE': [[75,200]],
                        }
            primer3_primary_result = primer3.bindings.designPrimers(seq_args, global_args)
            primer3_result_table_dict = {}
            for i in range(primer3_primary_result["PRIMER_PAIR_NUM_RETURNED"]):
                primer_id = str(i)
                for key in primer3_primary_result:
                    if primer_id in key:
                        info_tag = key.replace("_" + primer_id, "")# 要将每个信息中的数字和下划线去掉
                        try:# 就是把不同的引物对结果归到一起
                            primer3_result_table_dict[info_tag]
                        except:
                            primer3_result_table_dict[info_tag] = []
                        finally:
                            if isinstance(primer3_primary_result[key],float):
                                primer3_primary_result[key]=float('%.2f' % primer3_primary_result[key])
                                primer3_result_table_dict[info_tag].append(primer3_primary_result[key])
                            else:
                                primer3_result_table_dict[info_tag].append(primer3_primary_result[key])
            #挑出引物序列
            for k,v in primer3_result_table_dict.items():
                dict_left_seqs={}  
                dict_right_seqs={}  
                dict_intnl_seqs={}
                dict_prdct_sizes={}
                dict_compl_anyths={}
                dict_pair_penaltys={}  
                dict_right_RC_seqs={}       
                if k=="PRIMER_LEFT_SEQUENCE":
                    dict_left_seqs={"PRIMER_LEFT_SEQUENCE_"+str(count):v}
                    list_dict_left_seqs.append(dict_left_seqs)
                if k=="PRIMER_RIGHT_SEQUENCE":    
                    dict_right_seqs={"PRIMER_RIGHT_SEQUENCE_"+str(count):v} 
                    list_dict_right_seqs.append(dict_right_seqs)
                    vv_r_c=[]
                    for vv in v:
                        vv_r_c.append(Seq(vv,IUPAC.unambiguous_dna).reverse_complement())
                        dict_right_RC_seqs={"PRIMER_RIGHT_SEQUENCE_RC_"+str(count):vv_r_c}
                    list_dict_right_seqs_rc.append(dict_right_RC_seqs)                      
                if k=="PRIMER_INTERNAL_SEQUENCE":
                    if req.POST['PRIMER_PICK_INTERNAL_OLIGO'] == "0":
                        dict_intnl_seqs={"PRIMER_INTERNAL_SEQUENCE_"+str(count):None}
                    else:
                        dict_intnl_seqs={"PRIMER_INTERNAL_SEQUENCE_"+str(count):v}
                    list_dict_intnl_seqs.append(dict_intnl_seqs)
                if k=="PRIMER_PAIR_PRODUCT_SIZE":
                    dict_prdct_sizes={"PRIMER_PAIR_PRODUCT_SIZE_"+str(count):v}
                    list_dict_prdct_sizes.append(dict_prdct_sizes)
                if k=="PRIMER_PAIR_COMPL_ANY_TH":    
                    dict_compl_anyths={"PRIMER_PAIR_COMPL_ANY_TH_"+str(count):v} 
                    list_dict_compl_anyths.append(dict_compl_anyths)
                if k=="PRIMER_PAIR_PENALTY":
                    dict_pair_penaltys={"PRIMER_PAIR_PENALTY_"+str(count):v}
                    list_dict_pair_penaltys.append(dict_pair_penaltys)
            count+=1
        table_nums=len(list_dict_right_seqs)
        return render_to_response('multiprimer/multiseqview.html',{
                                    'list_dict_left_seqs':list_dict_left_seqs,
                                    'list_dict_right_seqs':list_dict_right_seqs,
                                    'list_dict_intnl_seqs':list_dict_intnl_seqs,
                                    'list_dict_prdct_sizes':list_dict_prdct_sizes,
                                    'list_dict_compl_anyths':list_dict_compl_anyths,
                                    'list_dict_pair_penaltys':list_dict_pair_penaltys,
                                    'list_dict_right_seqs_rc':list_dict_right_seqs_rc,
                                    'table_nums':table_nums,
                                    'inputseqs':inputseqs,
                                    'local':localip,
                                    'thisyear':thisyear
        },context_instance=RequestContext(req))