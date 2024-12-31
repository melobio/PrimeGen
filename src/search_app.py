import requests
import json
import os
# NER request function
def get_ner_entities(base_url, text):
    url = f"{base_url}/ner"
    payload = {"text": text}
    response = requests.post(url, json=payload)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"[Warning] NER request failed with status code {response.status_code}: {response.text} ")
        return False

# Embedding Search request function
def embedding_search(base_url, text, type_, th_score=0.6, top_k=10):
    url = f"{base_url}/embedding_search"
    payload = {
        "text": text,
        "type_": type_,
        "th_score": th_score,
        "top_k": top_k
    }
    response = requests.post(url, json=payload)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Embedding Search request failed with status code {response.status_code}: {response.text}")


# Synonym Check request function
def check_synonym(base_url, text_list):
    url = f"{base_url}/check_synonym"
    payload = {"text": text_list}
    response = requests.post(url, json=payload)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Check Synonym request failed with status code {response.status_code}: {response.text}")



# Get Meta Info request function
def get_meta_info(base_url, query_list):
    url = f"{base_url}/get_meta_info"
    payload = {"text": query_list}
    response = requests.post(url, json=payload)
    if response.status_code == 200:
        return response.json()
    else:
        raise Exception(f"Get Meta Info request failed with status code {response.status_code}: {response.text}")

def get_embedding_recommended(text, type, species_class=""):
    # Define request data
    data = {
            "text": text,
            "mode": type,
            "species_class": species_class,
            "top_k": 50,
            "top_r": 3,
            "score_threshold": 0.2
            }
    
    print("Sending data:", json.dumps(data, indent=4))
    
    SPECIES_SEARCH_HOST = os.getenv('SPECIES_SEARCH_HOST')
    SPECIES_SEARCH_PORT = os.getenv('SPECIES_SEARCH_PORT')
    # Send POST request
    response = requests.post(
                url=f'http://{SPECIES_SEARCH_HOST}:{SPECIES_SEARCH_PORT}/search/',
                headers={'Content-Type': 'application/json', 'accept': 'application/json'},
                json=data
                )
    # Check response status code
    if response.status_code == 200:
        # Print response content
        if type == "determine_species":
            print(response.json())
            print('+++++++++++++++')
            search_results = response.json()["species_results"]
            print("recommended search results sampling: ",search_results)
            return search_results

        elif type == "determine_cds":
            search_results = response.json()["cds_results"]
            return search_results

        elif type == "get_isolate2id":
            search_results = response.json()["isolate2id_results"]
            return search_results

    else:
        print(f"Request failed with status code {response.status_code}")