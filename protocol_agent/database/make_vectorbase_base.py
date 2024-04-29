import re
import json
from collections import defaultdict
import pickle
import numpy as np
import glob
from tqdm import tqdm
from sklearn.feature_extraction.text import TfidfVectorizer
from FlagEmbedding import BGEM3FlagModel

# text2vector
def text_to_vector_space(model,sentences):
    vectory_library = model.encode(sentences,batch_size=128, return_dense=True, return_sparse=True, return_colbert_vecs=False)
    return vectory_library


def chunk_text(text, chunk_size=2048):
    # Split the text by table tags
    segments = re.split(r'(<tablebegin>.*?<tableend>)', text, flags=re.DOTALL)
    chunks = []
    current_chunk = ''
    for segment in segments:
        if segment.startswith('<tablebegin>'):
            # Remove table tags
            table_text = segment.replace('<tablebegin>', '').replace('<tableend>', '')

            # If there's text before the table, append it as a chunk
            if current_chunk:
                chunks.append(current_chunk.strip())
                current_chunk = ''

            # Chunk the table text if it's too long
            for i in range(0, len(table_text), chunk_size):
                table_chunk = table_text[i:i+chunk_size].strip()
                if table_chunk:
                    chunks.append(table_chunk)
        else:
            # Process non-table text
            for char in segment:
                if len(current_chunk + char) < chunk_size:
                    current_chunk += char
                else:
                    chunks.append(current_chunk.strip())
                    current_chunk = char
    # Add any remaining text
    if current_chunk:
        chunks.append(current_chunk.strip())
    return chunks


#
def main(data_dict,chunk_length):
    # 创建向量库和映射
    file_library=[]
    txt_library=[]
    #构建embedding模型的向量数据
    for name,info_ in data_dict.items():
        text = info_["content"]
        chunks = [name]+[re.sub(r'-{2,}','-',chunk) for chunk in chunk_text(text, chunk_length) if chunk !=""]
        
        num=len(chunks)
        file_library+=[name]*num
        txt_library+=chunks

    return file_library,txt_library
    

if __name__=="__main__":
    model = BGEM3FlagModel("/home/YJhou/compspace/5_LLMWeights/Embedding_model/bge-m3")

    # text chunk
    chunk_length = 512  # 可以根据需要进行调整
    
    protocl_path = glob.glob("./protocol/*")
    API_path = glob.glob("./Opentrons_API_documents/*")

    file_path = API_path+protocl_path
    content_json = {}
    for file_ in tqdm(file_path):
        with open(file_, 'r') as f:
            content = f.read()
            content_json[file_.split('/')[-1].replace('.txt','')] = {}
            content_json[file_.split('/')[-1].replace('.txt','')]["content"] = content
        f.close()

    file_library,txt_library = main(content_json,chunk_length)    
    print('分块完成，开始编码')
    
    vector_library = text_to_vector_space(model,txt_library)
    with open("m3_vector_library.pickle", "wb") as f:
        pickle.dump(vector_library, f)
    f.close()
    np.save('m3_file_library.npy',np.array(file_library))
    np.save('m3_txt_library.npy',np.array(txt_library))
    
