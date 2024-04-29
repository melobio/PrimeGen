import openai
import re
import json
import os
import glob
import shutil
import Levenshtein
import time
import numpy as np
import pickle
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

app = FastAPI()

class Retrieval_():
    def __init__(self, embedding_path, reranking_path, txt_library_path, file_library_path, vector_library_path):
        self.LoadEmbeddingModel(embedding_path)
        self.LoadRerankingModel(reranking_path)
        self.LoadingVectorBase(txt_library_path, file_library_path, vector_library_path)

    
    def LoadEmbeddingModel(self, embedding_path):
        #now only support bge-m3
        from FlagEmbedding import BGEM3FlagModel
        self.EM_model = BGEM3FlagModel(embedding_path)

    def LoadRerankingModel(self, reranking_path):
        #now only support bge-m3
        from FlagEmbedding import FlagReranker
        self.reranker_model = FlagReranker(reranking_path)
    
    def LoadingVectorBase(self, txt_library_path, file_library_path, vector_library_path):
        #now only support pickle and npy files
        with open(vector_library_path,"rb") as f:
            self.doc_embeddings = pickle.load(f)
        f.close()
        self.file_library = np.load(file_library_path, allow_pickle=True)
        self.txt_library = np.load(txt_library_path, allow_pickle=True)

        

    def GetEmbedding(self, query):
        #now only support bge-m3
        self.embedding = self.EM_model.encode(query,
                                                return_dense=True,
                                                return_sparse=True,
                                                return_colbert_vecs=False)
        return self.embedding
    
    def ComputeEmbeddingSimilarityScore(self, query, top_k=10, score_threshold=0.5):
        similarity_score_list = []
        start_time = time.time()
        query_embeddings = self.GetEmbedding(query)


        for i in range(len(self.doc_embeddings["dense_vecs"])):
            dense_score = query_embeddings["dense_vecs"] @ self.doc_embeddings["dense_vecs"][i].T
            lexical_score = self.EM_model.compute_lexical_matching_score(query_embeddings['lexical_weights'],
                                                                    self.doc_embeddings['lexical_weights'][i])

            similarity_score = dense_score*0.8 + lexical_score*0.2 #+ multi_vector_score*0.4
            similarity_score_list.append(similarity_score)
        sorted_lst_with_idx = sorted(enumerate(similarity_score_list), key=lambda x: x[1], reverse=True)
        top_rank_score_and_idxs = sorted_lst_with_idx[:top_k]
        self.embedding_results = []
        for i,score in top_rank_score_and_idxs:
            if score > score_threshold:
                self.embedding_results.append({'text':self.txt_library[i], 'file':self.file_library[i], 'score':score,'index':i})
                                     
        print("embedding search time:", time.time()-start_time)
        return self.embedding_results
    
    def ComputeRerankScore(self, query,top_r=2):
        retrieval_pairs = [[query,result['file']+'\n'+result['text']] for result in self.embedding_results]
        scores = self.reranker_model.compute_score(retrieval_pairs, normalize=True)
        scored_results = [{'text': result['text'], 'file': result['file'], 'score': score} for result, score in zip(self.embedding_results, scores)]
        self.reranked_results = sorted(scored_results, key=lambda x: x['score'], reverse=True)[:top_r]
        return self.reranked_results
    
    def LLMSearch(self,):
        #GPT4 socre is coming soon
        pass

    def retrieval_pipeline(self, query, top_k=10, score_threshold=0.5, top_r=2):
        embedding_results = self.ComputeEmbeddingSimilarityScore(query, top_k, score_threshold)
        reranked_results = self.ComputeRerankScore(embedding_results, top_r)
        #llmscore_results
        return reranked_results
    
class Query(BaseModel):
    text: str
    top_k: str
    top_r: str
    score_threshold: str


# 全局变量初始化
retrieval = None

@app.on_event("startup")
async def startup_event():
    global retrieval
    embedding_path = "/home/YJhou/compspace/5_LLMWeights/Embedding_model/bge-m3"
    reranking_path = "/home/YJhou/compspace/5_LLMWeights/Embedding_model/bge-reranker-v2-m3"
    txt_library_path = ""
    file_library_path = ""
    vector_library_path = ""
    retrieval = Retrieval_(embedding_path, reranking_path, txt_library_path, file_library_path, vector_library_path)



@app.post("/search/")
def search(query: Query):
    results = retrieval.retrieval_pipeline(query.text, int(query.top_k), int(query.top_r), float(query.score_threshold))
    return {"retrieval_results":results}

if __name__=="__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8008)

    """#Use Cases
    curl -X 'POST' \
    'http://127.0.0.1:8000/search/' \
    -H 'accept: application/json' \
    -H 'Content-Type: application/json' \
    -d '{
    "text": "example query",
    "top_k": 10,
    "top_r": 2,
    "score_threshold": 0.5
    }'
    """
    
