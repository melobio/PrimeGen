import os 
import json
from pathlib import Path
from typing import Annotated, Union
from tqdm import tqdm
import typer
from swift.llm import (
     get_template,  ModelType,
    get_default_template_type,  inference_vllm, get_vllm_engine
)


app = typer.Typer(pretty_exceptions_show_locals=False)

def load_model_and_template(
        model_id_or_path: Union[str, Path]
):
    model_type = ModelType.phi3_vision_128k_instruct
    llm_engine = get_vllm_engine(model_type, model_id_or_path = model_id_or_path, tensor_parallel_size=4)
    template_type = get_default_template_type(model_type)
    template = get_template(template_type, llm_engine.hf_tokenizer)
    llm_engine.generation_config.max_new_tokens = 256
    return llm_engine, template

def run(model, template,img):
    images = [img] 
    request_list = [{'query': 'Is the OT-2 Gen-2 robotic arm gripper holding a pipette tip, and if so, what is the condition of the pipette tip? Also, is the heated lid of the Thermocycler Module open or closed, and how many pipette tips are held by the robotic arm? In the current image, is the trash can empty or full? [Holding/Released, Undamaged/Damaged/None, Open/Closed, One/Eight/None, Empty/Full]', 'images': images}]
    resp_list = inference_vllm(model, template, request_list)
    resp = resp_list[0]

    return resp['response']

    

def finding_n(String,Substr,times): 
    '''
    Find the index of the times-th occurrence of Substr in a String.
    '''
    nemo = 0
    for i in range(1,times+1):
        nemo = String.find(Substr,nemo) + 1
        if nemo == 0:
            return -1
    return nemo-1
    
def save_json(val_metrics,json_path="val_metrics.json"):
    data = json.dumps(val_metrics)
    f2 = open(json_path, 'w')
    f2.write(data)
    f2.close()

def load_json(json_path="val_metrics.json"):
    f = open(json_path, 'r')
    content = f.read()
    res = json.loads(content)
    f.close()
    return res
class Predict:    
    def __init__(self, model, template):
        self.model = model
        self.template = template
    
    def pickUpTip(self,action_dir,question):
        
        preds = []
        total_img = len(os.listdir(action_dir))        
        for i in tqdm(range(0, total_img), total=total_img):
            name = sorted(os.listdir(action_dir))[i]
            
            j = os.path.join(action_dir,name)
            one_result = run(self.model,self.template,j)

            # For the issue of pipette tip presence, it is 'A' if present, otherwise 'B'
            if  question == 'single':
                print(one_result)
                if ("Holding" in one_result.split(',')[0]):
                    preds.append('A')
                else:
                    preds.append('B')
                print(preds)
            # or the issue of pipette tip presence, it is 'A' if present, otherwise 'B'      
            elif question == 'multi':
                if 'Eight' in one_result.split(',')[3]:
                    preds.append('A')
                else:
                    preds.append('B')
                print(preds)

        str_pred = "".join(map(str, preds))  
        idx = finding_n(str_pred,"AAA",1)
        print(str_pred)
        print(idx)
        if len(preds) == 0:
            return False
        elif idx == 0 or idx == -1:
            return False
        else:
            return True
            
    def dropTip(self,action_dir,question):
        preds = []
        total_img = len(os.listdir(action_dir))
        for i in tqdm(range(0, total_img), total=total_img):
            name = sorted(os.listdir(action_dir))[i]
            j = os.path.join(action_dir,name)
            one_result = run(self.model,self.template,j)
            print(one_result)
            # For the issue of pipette tip presence, it is 'B' if present, otherwise 'a'
            if  question == 'single':
                if "Released" in one_result.split(',')[0]: #Not Holding
                    preds.append('B')
                else:
                    preds.append('A')
                print(preds)

            elif question == 'multi':
                if 'Eight' in one_result.split(',')[3]:
                    preds.append('A')
                else:
                    preds.append('B')
                print(preds)

        str_pred = "".join(map(str, preds))
        idx = str_pred.rfind("BBB")
        print(str_pred)
        print(idx)
        if len(preds) == 0:
            return False
        elif idx == len(preds) - 1:
            return False
        elif idx == -1:
            return True
        else:
            return True
    
    def PCR_Open(self,action_dir):
        question = 'lid'
        preds = []
        total_img = len(os.listdir(action_dir))
        for i in tqdm(range(0, total_img), total=total_img):
            name = sorted(os.listdir(action_dir))[i]
            j = os.path.join(action_dir,name)
            one_result = run(self.model,self.template,j)
            print(one_result)
            if 'Open' in one_result.split(',')[2]:
                preds.append('B')
            else:
                preds.append('A')
            print(preds)
        
        str_pred = "".join(map(str, preds))
        idx = str_pred.rfind("BBB") 

        if len(preds) == 0:
            return False

        elif idx == len(preds) - 1: 
            return False
        elif idx == -1:
            return True
        else:
            return True
        
    def PCR_Close(self,action_dir):
        question = 'lid'
        preds = []
        total_img = len(os.listdir(action_dir))

        for i in tqdm(range(0, total_img), total=total_img):
            name = sorted(os.listdir(action_dir))[i]
            j = os.path.join(action_dir,name)
            one_result = run(self.model,self.template,j)
            print(one_result)

            if 'Open' in one_result.split(',')[2]:
                preds.append('B')
            else:
                preds.append('A')
            print(preds)
        str_pred = "".join(map(str, preds))  
        print(f"string{str_pred}")
        idx = finding_n(str_pred,"AAA",1)
        print(str_pred)
        print(idx)
        if len(preds) == 0:
            return False
        elif idx == 0 or idx == -1:
            return False
        else:
            return True

@app.command()
def main(
        model_dir: Annotated[str, typer.Argument(help='Path to model')],
        save_file: Annotated[str, typer.Argument(help='Path to save jsonl')]
):
    
    model, template = load_model_and_template(model_dir)             
    #The user needs to provide the image addresses for consecutive frames.
    root_dir = ""
    actions = ["close_lid","open_lid","drop","drop_8","pickup_8","pickup"]

    acc = {}
    multi_ins = Predict(model, template)

###############################################################
    
    
    for act in actions:
        action_dir = os.path.join(root_dir,act)
        all = 0
        pred = 0
        
        
        if act == "close_lid":
            print(f"Processing : {act}")
            for p in os.listdir(action_dir): # [1 2 3 4 5 6 7]
                if p=='2' or p=='8':
                    continue
                print(action_dir+'/'+p)
                all += 1
                res = multi_ins.PCR_Close(os.path.join(action_dir,p))
                pred+=int(res)
                
            print(f"【close_lid】: total = {all} , correct = {pred}  acc = {pred/all}")
       
        #if act == "open_lid":
        elif act == "open_lid":
            print(f"Processing : {act}")
            for p in sorted(os.listdir(action_dir)):
                if p=='2' or p=='8':
                    continue
                print(action_dir+'/'+p)
                all += 1
                res = multi_ins.PCR_Open(os.path.join(action_dir,p))
                pred+=int(res)
            print(f"【open_lid】: total = {all} , correct = {pred}  acc = {pred/all}")
       
        elif act == "pickup":
            print(f"Processing : {act}")
        #if act == "pickup":
            for p in os.listdir(action_dir):
                if p=='2' or p=='4' or p=='6' or p=='8':
                    continue
                print(action_dir)
                print(p)
                #if p=='1' or p=='2':
                #    continue
                all += 1
                res = multi_ins.pickUpTip(os.path.join(action_dir,p),'single')
                pred+=int(res)
            print(f"【pickup】: total = {all} , correct = {pred}  acc = {pred/all}")
        
        elif act == "pickup_8":
            print(f"Processing : {act}")
            for p in os.listdir(action_dir):
                if p=='2' or p=='4' or p=='6' or p=='8':
                    continue
                all += 1
                res = multi_ins.pickUpTip(os.path.join(action_dir,p),'multi')
                pred+=int(res)
            print(f"【pickup_8】: total = {all} , correct = {pred}  acc = {pred/all}")
        elif act == "drop":
            print(f"Processing : {act}")
            for p in os.listdir(action_dir):
                if p=='2' or p=='8':
                    continue
                all += 1
                res = multi_ins.dropTip(os.path.join(action_dir,p),'single')
                pred+=int(res)
            print(f"【drop】: total = {all} , correct = {pred}  acc = {pred/all}")
        elif act == "drop_8":
            print(f"Processing : {act}")
            for p in os.listdir(action_dir):
                if p=='2' or p=='8':
                    continue
                all += 1
                res = multi_ins.dropTip(os.path.join(action_dir,p),'multi')
                pred+=int(res)
            print(f"【drop_8】: total = {all} , correct = {pred}  acc = {pred/all}")
        else:
            continue
        tmp = acc.get(act,[])
        print(act)
        print(all)
        print(pred)
        tmp.append(all)
        tmp.append(pred)
        tmp.append(pred/all)
        acc[act] = tmp

    name_save = 'Real_Phi3_acc.json'
    save_json(acc,json_path=os.path.join(save_file,name_save))


if __name__ == '__main__':
    app()
