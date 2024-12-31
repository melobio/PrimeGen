from typing import Annotated
import typer
import json
from openai import OpenAI

import base64
from tqdm import tqdm
import csv
import nltk
#nltk.download('wordnet')
from nltk.translate.bleu_score import sentence_bleu, SmoothingFunction
from rouge import Rouge
from openai import AzureOpenAI

from infer_prompt import score_system_prompt, user_prompt_template

rouge = Rouge()
app = typer.Typer(pretty_exceptions_show_locals=False)
gpt4_client = AzureOpenAI(
                azure_endpoint = "***",
                api_key="***",
                api_version="***"
                )

# use base64
def get_image_base64(img_path):
    with open(img_path, 'rb') as f:
        img_base64 = base64.b64encode(f.read()).decode('utf-8')
    image_url = f'data:image/jpeg;base64,{img_base64}'
    return image_url


def get_response(query, img_path, model_type, client):
    image_url = get_image_base64(img_path)
    messages = [{
            'role': 'user',
            'content': [
                {'type': 'image_url', 'image_url': {'url': image_url}},
                {'type': 'text', 'text': query.replace("<image>","")},
                ]
            }]

    resp = client.chat.completions.create(
                model=model_type,
                messages=messages,
                max_tokens=512,
                seed=42)

    response = resp.choices[0].message.content
    return response

def gpt4_score(ture_answer, prediction):
    instruction = user_prompt_template(ture_answer, prediction)
    messages = [
            { "role": "system", "content": score_system_prompt},
            { "role": "user", "content": [
                {   
                    "type": "text",
                    "text": instruction
                    },
                ] }
            ]
    response = gpt4_client.chat.completions.create(
                model="gpt-4o",
                response_format={"type": "json_object"},
                messages=messages,
                temperature=0.8,
                top_p=0.95,
                max_tokens=16384
                )
    output = (response.choices[0].message.content)
    description_response_dict = json.loads(output)
    score = description_response_dict['similarity_score']
    justification = description_response_dict['justification']
    return score, justification

def calculate_scores(true_answer, response):
    chencherry = SmoothingFunction()
    bleu_score = sentence_bleu([true_answer.split()], response.split(), smoothing_function=chencherry.method1)
    rouge_score = rouge.get_scores(response, true_answer)[0]
    rouge_l = rouge_score['rouge-l']['f']

    meteor='none'
    llm_score, llm_justification = gpt4_score(true_answer, response)
    return bleu_score, rouge_l, meteor, llm_score, llm_justification  

@app.command()
def main(
        data_path: Annotated[str, typer.Argument(help='Path to JSONL')],
        save_file: Annotated[str, typer.Argument(help='Path to save jsonl')]
):  
    client = OpenAI(
                    api_key='EMPTY',
                    base_url='http://localhost:8181/v1',
                    )

    model_type = client.models.list().data[0].id
    print(f'model_type: {model_type}')

    with open(save_file, mode='w', newline='', encoding='utf-8') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['query', 'img_path', 'true_answer', 'prediction', 'BLEU', 'ROUGE-L', 'METEOR','LLMscore','LLMjustification'])
        with open(data_path, 'r', encoding='utf-8') as f:
            lines = list(f)
            for line in tqdm(lines, total=len(lines), desc="Processing"):
                conversation = json.loads(line)

                query = conversation['query']
                img_path = conversation['images'][0]
                true_answer = conversation['response']
                print(query)
                print("【True】:", true_answer)
                response = get_response(query, img_path, model_type, client)
                print("【Prediction】:",response)
                bleu_score, rouge_l, meteor, llm_score, llm_justification = calculate_scores(true_answer, response)
                csv_writer.writerow([query, img_path, true_answer, response, bleu_score, rouge_l, meteor, llm_score, llm_justification])

                print(f'【BLEU Score】: {bleu_score}')
                print(f'【ROUGE-L Score】: {rouge_l}')
                print(f'【LLM Score】: {llm_score}.【LLM Justification】: {llm_justification}')
                print('+++++++'*6)


if __name__ == '__main__':
    app()