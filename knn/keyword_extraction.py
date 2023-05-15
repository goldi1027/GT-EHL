import string
import re
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer
from nltk.stem import WordNetLemmatizer




with open('original.txt') as f:
    lines = f.readlines()
f.close()

for i in range (len(lines)):
    # convert to lower case
    lines[i] = lines[i].lower()
    # remove numbers
    lines[i] = re.sub(r'[0-9]+', '', lines[i])
    # remove remove punctuation
    lines[i] = lines[i].translate(str.maketrans('', '', string.punctuation))
    # remove whitespaces
    lines[i] = lines[i].strip()
    # tokenization
    lines[i] = word_tokenize(lines[i])
    # remove stop words
    stop_words = set(stopwords.words("english"))
    lines[i] = [j for j in lines[i] if not j in stop_words]
    # word stemming
    stemmer= PorterStemmer()
    lines[i] = [stemmer.stem(j) for j in lines[i]]
    # word lemmatization
    lemmatizer=WordNetLemmatizer()
    lines[i] = [lemmatizer.lemmatize(j) for j in lines[i]]
    # merge back to a single string using space to separate key words
    lines[i] = ' '.join(lines[i])

output = '\n'.join(lines)

file = open("cleaned.txt", "w")
file.write(output)
file.close()

print(output)