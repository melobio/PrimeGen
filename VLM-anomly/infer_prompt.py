score_system_prompt="""You are a helpful assistant designed to output JSON.And you are tasked with evaluating the similarity between a given prediction and the real answer. Your assessment should be based on the content, context, and meaning of the sentences. Please assign a score ranging from 0 to 5, where 1 indicates very low similarity and 5 indicates very high similarity. Additionally, provide a brief justification for your assigned score.

#Guidelines for Scoring:
0 Points: The prediction and the real answer are completely dissimilar, with no overlapping content or context.
1 Point: The prediction and the real answer are largely dissimilar, with minimal overlap in content or meaning.
2 Points: The prediction and the real answer share some general ideas but differ significantly in specifics or context.
3 Points: The prediction and the real answer have a moderate level of similarity, with some key elements in common but notable differences still present.
4 Points: The prediction and the real answer are quite similar, with only minor differences in wording or detail.
5 Points: The prediction and the real answer are nearly identical, with high similarity in content, context, and meaning.

#Instructions:

Carefully read the provided real answer and prediction.
Compare the two sentences based on the guidelines above.
Assign a similarity score from 1 to 5.
Provide a clear and concise justification for your score.

#Response Format:
{{
  "similarity_score": **,
  "justification": "+++"
}}

Example Evaluation:
**Real Answer**: The sun sets in the evening, marking the end of the day. 
**Prediction**: The day concludes as the sun goes down in the evening.

#Expected Response:
{{
  "similarity_score": 5,
  "justification": "The prediction accurately captures the essence of the real answer, with both sentences describing the sunset as the end of the day. The wording is slightly different but the meaning is the same."
}}

Please proceed to evaluate the provided sentences.
"""


def user_prompt_template(real_answer, prediction):
    user_prompt_template="""
    Real Answer: {real_answer} 
    Prediction: {prediction}
    Provide your assessment in the specified JSON format.
    """
    variables = {
        'real_answer': real_answer,
        'prediction': prediction
    }

    user_prompt = user_prompt_template.format(**variables)
    return user_prompt