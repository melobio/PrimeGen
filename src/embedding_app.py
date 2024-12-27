import time
import numpy as np
import pickle
# from fastapi import FastAPI, HTTPException, Depends
from pydantic import BaseModel
import pandas as pd
import torch
import argparse
# app = FastAPI()
from FlagEmbedding import BGEM3FlagModel
import openai
import json
import traceback
import os

from llm_utils.utils import get_llm_chat_completion


class Retrieval_():
    def __init__(self, embedding_path, reranking_path, txt_library_path, vector_library_path):
        self.LoadEmbeddingModel(embedding_path)
        # self.LoadRerankingModel(reranking_path)
        self.LoadingVectorBase(txt_library_path, vector_library_path)

    # 加载Embedding模型
    def LoadEmbeddingModel(self, embedding_path):
        #暂时只支持bge-m3模型
#         from FlagEmbedding import BGEM3FlagModel
        self.EM_model = BGEM3FlagModel(embedding_path)
    
    # 加载Rerank模型
    # def LoadRerankingModel(self, reranking_path):
    #     #暂时只支持bge-m3模型
    #     from FlagEmbedding import FlagReranker
    #     self.reranker_model = FlagReranker(reranking_path)
    
    # 定义加载向量和文本
    def LoadingVectorBase(self, txt_library_path, vector_library_path):
        #暂时pickle和npy文件
        with open(vector_library_path,"rb") as f:
            self.doc_embeddings = pickle.load(f)
        f.close()
        #self.file_library = np.load(file_library_path, allow_pickle=True)
#         self.txt_library = np.load(txt_library_path, allow_pickle=True)
        des_df = pd.read_csv(txt_library_path)
        self.txt_library =list(des_df['description'])

    # 定义嵌入方法
    def GetEmbedding(self, query):
        #暂时只支持bge-m3模型
        self.embedding = self.EM_model.encode(query,
                                                return_dense=True,
                                                return_sparse=False,
                                                return_colbert_vecs=False)
        return self.embedding
    
    # 定义Embedding检索方法——用pytorch加速矩阵运算
    def ComputeEmbeddingSimilarityScore(self, query, top_k=50, score_threshold=0.5):
        similarity_score_list = []
        start_time = time.time()
        query_embeddings = self.GetEmbedding(query)

        similarity_score_list = torch.matmul(torch.tensor(query_embeddings["dense_vecs"]).cuda(), torch.tensor(self.doc_embeddings["dense_vecs"]).cuda().T)
        topk_scores, topk_indices = torch.topk(similarity_score_list, top_k, dim=0, largest=True, sorted=True)
        scores_list = topk_scores.cpu().tolist()
        index_list = topk_indices.cpu().tolist()
        self.embedding_results = []
        for i,score in zip(index_list,scores_list):
            if score > score_threshold:
                self.embedding_results.append({'text':self.txt_library[i], 'score':score,'index':i})
        return self.embedding_results

    # 定义Rerank方法
    def ComputeRerankScore(self, query, top_r=10):
        retrieval_pairs = [[query,result['text']] for result in self.embedding_results]
        scores = self.reranker_model.compute_score(retrieval_pairs, normalize=True)
        scored_results = [{'text': result['text'], 'score': score} for result, score in zip(self.embedding_results, scores)]
        self.reranked_results = sorted(scored_results, key=lambda x: x['score'], reverse=True)[:top_r]
        return self.reranked_results

    def retrieval_pipeline(self, query, top_k=50, top_r=10, score_threshold=0.5):
        #先进行Embedding检索
        start_time=time.time()
        embedding_results = self.ComputeEmbeddingSimilarityScore(query, top_k, score_threshold)
#         print('embedding_results',embedding_results)
        print("Embedding耗时：",time.time()-start_time)
        
        return embedding_results
#         #再进行Rerank排序
#         start_time=time.time()
#         reranked_results = self.ComputeRerankScore(query, top_r)
#         print('reranked_results',reranked_results)
#         print("Reranking耗时：",time.time()-start_time)
#         return reranked_results
    def retrieval_pipeline2(self, query, top_k=50, top_r=10, score_threshold=0.5):
        #先进行Embedding检索
        start_time=time.time()
        embedding_results = self.ComputeEmbeddingSimilarityScore(query, top_k, score_threshold)
#         print('embedding_results',embedding_results)
        print("Embedding耗时：",time.time()-start_time)
#         #再进行Rerank排序
        start_time=time.time()
        reranked_results = self.ComputeRerankScore(query, top_r)
        print('reranked_results',reranked_results)
        print("Reranking耗时：",time.time()-start_time)
        return reranked_results

#定义入参解析格式
class Query(BaseModel):
    text: str
    top_k: int
    top_r: int
    score_threshold: float

# # 全局变量存储模型和数据
# retrieval_instance = None
# # 启动服务先加载内容
# @app.on_event("startup")
# def load_models():
#     global retrieval_instance
#     # 预设定Embedding模型
#     embedding_path = "./bge-m3"
#     # 预设定Rerank模型
#     reranking_path = "./bge-reranker-v2-m3/"
#     # 文本地址
#     txt_library_path = "./code_block_description.csv"#alil_txt_library.npy"
#     # 向量地址
#     vector_library_path = "./data/all_m3_vector_library.pickle"#all_m3_vector_library.pickle"
    
#     # 构建检索方法
#     retrieval_instance = Retrieval_(embedding_path, reranking_path, txt_library_path, vector_library_path)

if __name__=="__main__":
    openai.api_type = "azure"
    openai.api_base = "https://xmgi-chat8.openai.azure.com/"
    openai.api_version = "2024-02-15-preview"
    openai.api_key = "2d4b3509a9904372be2a75f45a12203f"
    
    parser = argparse.ArgumentParser(description='text embedding ')
    parser.add_argument("--query_file_path", type=str, help="input file")
    args = parser.parse_args()
    query_file_path = args.query_file_path
    
    
    # 预设定Embedding模型
    embedding_path = "/reference_data/bge-m3"
    # 预设定Rerank模型
    reranking_path = "/reference_data/bge-reranker-v2-m3/"
    # 文本地址
    txt_library_path_1 = "/reference_data/block_description_filter.csv"#alil_txt_library.npy"
    txt_library_path_2 = "/reference_data/block_name_description.csv"
    # 向量地址
    vector_library_path_1 = "/reference_data/all_code_block_vector_library.pickle"#all_m3_vector_library.pickle"
    vector_library_path_2 = "/reference_data/all_name_vector_library.pickle"#all_m3_vector_library.pickle"
    # 构建检索方法
     
    retrieval_2 = Retrieval_(embedding_path, reranking_path, txt_library_path_2, vector_library_path_2)
    #通过步骤向量和描述向量两种方式进行多路召回
    query_df = pd.read_csv(query_file_path)
    df_all = pd.read_csv(txt_library_path_1)
    
    name_list = list(query_df['name'])
    text_list = list(query_df['text'])
    
    code_name_list_1 = []
    temp_sorted_list = []
    dict_list = []
    for nn in name_list:
        if '(' in nn:
            nn = nn.split('(')[0]
        temp_res = retrieval_2.retrieval_pipeline(nn, 10, 10, 0.1)
        print('nn',nn)
        print('temp_res',temp_res)
        temp_list =[]
        temp_n =[]
        temp_dict_list = []
        for xx in temp_res:
            temp_dict = {}
            temp_n.append(xx['text'])
            ind =list(df_all['new_name']).index(xx['text'])
            temp_list.append(list(df_all['description'])[ind])
            temp_dict['名称'] = xx['text']
            temp_dict['描述'] = list(df_all['description'])[ind]
            temp_dict_list.append(temp_dict)
        dict_list.append(temp_dict_list)
        code_name_list_1.append(temp_list)
        temp_sorted_list.append(temp_n)
    find_name_list = []
    for ii in range(len(name_list)):
        query_dict = {}
        query_dict['名称'] = name_list[ii]
        query_dict['描述'] = text_list[ii]
        target_dict_list = dict_list[ii]
        txt_file = 'query字典为：'+str(query_dict)+'  target字典库为: '+', '.join([str(d) for d in target_dict_list])
        print('query字典',query_dict)
        print('target字典库',target_dict_list)

        #通过GPT4来选择描述最相近的代码块
        experiment_type_prompt = (f"""
        用户会给你一个target字典和一个query字典库，字典中都有"名称"和"描述"两个关键字，你需要在target字典库中找出名称与query名称是同类型的并且描述与query的描述语义最相近的一个字典，并以JSON进行返回。
        名称是同类型的判断标准为："连接产物纯化"与"连接产物纯化2"是同类型的，另外打断后磁珠双选和打断产物纯化是同类型,多重扩增产物纯化和PCR产物纯化是同类型，打断末端修复与酶切打断是同类型。
        描述最相近的的判断标准为：
        1. 所配置的反应液要是同一种；
        2. 移液的次数和所使用到的溶液数量要相同，DNA Clean Beads与En-Beads是同一种溶液，En-TE与TE Buffer是同一种溶液；
        3. 当名称是接头连接反应时，如果在他的描述中，吸取接头连接反应液之后没有再加入TE_buffer ，则他是与接头连接反应2相同； 
        4.当名称是连接产物纯化时,如果在他的描述中，出现了两次加入160μL 80%乙醇则他是与连接产物纯化2相同；
       
        以下是一个具体的例子：
        query字典为 {{"名称": "接头连接反应", "描述": "在冰上配制接头连接反应混合液：24μL Ligation Buffer Mix和1μL Ligation Enzyme。在末端修复产物的PCR管中加入5μL对应的Adapters Mix，涡旋震荡3次，每次3秒，瞬时离心。用移液器缓慢吸取25μL配制好的接头连接反应混合液加入样本管中，涡旋震荡6次，每次3秒，瞬时离心。将PCR管置于PCR仪上，按照以下条件进行反应：23℃ 20 min，4℃ Hold。"}}
        target字典库 [{{
            "名称": "接头连接反应1",
            "描述": "开始时，在冷环境下准备接头连接反应液，混合Ligation Buffer 23.4 μL与DNA Ligase 1.6 μL，总体积为25 μL。向PCR管添加5 μL的MGIEasy DNA Adapters后，进行三次快速涡旋，每次持续3秒，接着进行快速离心，以集中管底的液体。随后，慢慢吸取已配好的25 μL接头连接反应液到PCR管中，进行六次再涡旋，每次也是3秒，再次短暂离心确保液体集中。放入PCR仪进行预设反应，完成后短暂离心收集反应液。最终，加入20 μL TE Buffer至PCR管中，调整总体积至100 μL，并转移至新的1.5 mL离心管完成最后步骤。"
        }}, {{'名称': '接头连接反应2', '描述': '首先，将UDB Adapters和Fast Ligation Buffer从储存中取出至常温，使其解冻后充分涡旋。彻底摇动Ad Ligase 5-10次以混合均匀，瞬间离心后放置于低温环境。在冷环境中根据需求配制接头连接反应液：混合Fast Ligation Buffer 23μL、Ad Ligase 5μL和Ligation Enhancer 2μL至总体积30μL。先于PCR管中加入5μL UDB Adapters，随后缓慢注入30μL预配的反应液，六次涡旋，每次3秒，快速离心后集中管底液体。置于PCR仪进行设定的反应程序：在25℃下反应10分钟，之后保持在4℃。'}},{{'名称': '连接产物纯化3', '描述': '1. 混匀En-Beads，分别吸取20μL En-Beads和20μL En-TE至样本管中，用移液器轻轻吹打至少10次至所有磁珠悬浮。2. 室温孵育5 min。3.将PCR管瞬时离心，再置于磁力架上静置2～5 min至液体澄清，小心吸取上清并丢弃。4. 保持PCR管固定于磁力架上，加入160μL 80%乙醇漂洗磁珠及管壁，静置30s，小心吸取上清并丢弃。5. 重复步骤4一次，尽量吸干管内液体。6. 保持PCR管固定于磁力架上，打开管盖，室温干燥，直至磁珠表面无反光、无开裂。7. 将PCR管从磁力架上取下，加入20μL En-TE进行DNA洗脱，用移液器轻轻吹打至少10次至所有磁珠悬浮。8. 吸取20μL En-Beads至步骤7中的PCR管中，用移液器轻轻吹打至少10次至所有磁珠悬浮。9. 室温孵育5 min。10. 将PCR管瞬时离心，再置于磁力架上静置2～5 min至液体澄清，小心吸取上清并丢弃。11. 保持PCR管固定于磁力架上，加入160μL 80%乙醇漂洗磁珠及管壁，静置30s，小心吸取上清并丢弃。12. 重复步骤11一次。尽量吸干管内液体。13. 保持PCR管固定于磁力架上，打开管盖，室温干燥，直至磁珠表面无反光、无开裂。14.将PCR管从磁力架上取下，加入27μL En-TE进行DNA洗脱，用移液器轻轻吹打至少10次至所有磁珠悬浮。15. 室温孵育5 min。16. 将PCR管瞬时离心，再置于磁力架上静置2～5 min至液体澄清，用移液器小心吸取25μL上清液至新的0.2 mL PCR管。'}},{{'名称': '末端修复反应1', '描述': ' 取DNA样本至0.2mL PCR管中，不足40 μL部分用TE Buffer补足。在冰上配制末端修复反应液（ERAT Buffer 7.1 μL, ERAT Enzyme Mix 2.9 μL）。用移液器吸取10 μL配制好的末端修复反应液加入磁珠片段筛选或纯化后的40 μL样本中，涡旋震荡3次，每次3s，瞬时离心将反应液收集至管底。将PCR管置于PCR仪上，按照表9中的条件进行反应。反应结束后，瞬时离心将反应液收集至管底。'}}]

        返回的JSON为：{{
            "名称": "接头连接反应2",
            "描述": "首先，将UDB Adapters和Fast Ligation Buffer从储存中取出至常温，使其解冻后充分涡旋。彻底摇动Ad Ligase 5-10次以混合均匀，瞬间离心后放置于低温环境。在冷环境中根据需求配制接头连接反应液：混合Fast Ligation Buffer 23μL、Ad Ligase 5μL和Ligation Enhancer 2μL至总体积30μL。先于PCR管中加入5μL UDB Adapters，随后缓慢注入30μL预配的反应液，六次涡旋，每次3秒，快速离心后集中管底液体。置于PCR仪进行设定的反应程序：在25℃下反应10分钟，之后保持在4℃。"
        }}
        """)
        experiment_info_response = get_llm_chat_completion(
            messages=[
                {"role": "system", "content": experiment_type_prompt},
                {"role": "user", "content": txt_file}
            ],
            response_format={"type": "json_object"},
            temperature=0.2,
        )
        


        print('experiment_info_response', experiment_info_response.choices[0].message.content)
        # experiment_process = json.load(experiment_info_response['choices'][0]['message']['content'])
        experiment_process = experiment_info_response.choices[0].message.content
        if experiment_process:
            experiment_process = json.loads(experiment_process)
            find_name_list.append(experiment_process['名称'])


    max_name_list = [list(df_all['code_block'])[list(df_all['new_name']).index(na)]for na in find_name_list]    
    print('max_name_list',max_name_list)
    query_df.insert(1,'code_block',max_name_list)
    query_df.to_csv(query_file_path,index=False)
