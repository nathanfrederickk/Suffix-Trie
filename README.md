# Suffix-Trie
Uses a suffix trie to determine the longest common substring between two sentences. The result will be
a list containing: [longest common substring, the similarity score for submission1 expressed as percentage, 
                    the similarity score for submission2 expressed as percentage].
The percetage is rounded to the nearest integer

For example:
>>> submission1 = ’the quick brown fox jumped over the lazy dog’
>>> 
>>> submission2 = ’my lazy dog has eaten my homework’
>>> 
>>> compare_subs(submission1, submission2)
>>> 
>>> [’ lazy dog’, 20, 27]
