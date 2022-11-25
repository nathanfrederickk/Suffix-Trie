class Node:
    def __init__(self, type = 1) -> None:
        """
        Initialize a node to be used in a Trie. 

        :Input:
        type: an integer to determine where the string is coming from. 1 if it
              is coming from string 1, and 2 if it is coming from string 2

        :Output, return or postcondition: Returns a node with 28 spaces, 1 for
                                          the end marking, 1 for space, and 26 for
                                          lower case alphabets

        :Time complexity: O(1)
        :Aux space complexity: O(1)
        :Citation: FIT2004 2020sem02 Lecture11 Trie
        """
        self.link = [None] * 28
        self.type = type

class SuffixTrie:
    def __init__(self) -> None:
        """
        Initialize a suffix Trie with a root Node.

        :Output, return or postcondition: Returns an empty Trie with
                                          a root.

        :Time complexity: O(1)
        :Aux space complexity: O(1)
        :Citation: FIT2004 2020sem02 Lecture11 Trie
        """

        self.root = Node(0)
        self.path = []

    def _get_index(self, char):
        """
        Turns a char which is a one letter string to an index
        to be positioned into the link inside the node.

        :Input:
        char: a single letter string

        :Output, return or postcondition: Returns the index of where
                                          the string will be stored in the 
                                          link node.

        :Time complexity: O(1)
        :Aux space complexity: O(1)
        """
        return max(ord(char) - 97 + 2, 1)
    
    def insert(self, key):
        """
        Insert the first key and its suffixes to the Trie
        by using a recursion.

        :Input:
        key: a string which will be inserted to the Trie along
             with its suffixes.

        :Output, return or postcondition: The key is inserted to the string
                                          along with its suffixes.

        :Time complexity: O(N^2), where N is the length of key
        :Aux space complexity: O(N^2)
        :Citation: FIT2004 2020sem02 Lecture11 Trie
        """
        for i in range(len(key) + 1):
            current = self.root
            self.insert_aux(current, key, i)

    def insert_aux(self, current, key, i) -> None:
        """
        Recursion to insert the key[i::] into the Trie.

        :Input:
        current : a node object which is the self.root
        key     : the string which is inserted into the Trie
        i       : an integer to determine which index of the key will be inserted

        :Output, return or postcondition: The key[i::] is inserted to the string
                                          

        :Time complexity: O(N), where N is the length of key
        :Aux space complexity: O(N)
        :Citation: FIT2004 2020sem02 Lecture11 Trie
        """
        if i == len(key):
            if current.link[0] is not None:
                current = current.link[0]
            else:
                current.link[0] = Node()
                current = current.link[0]

            return

        else:
            index = self._get_index(key[i])

            if current.link[index] is not None:
                current = current.link[index]
            else:
                current.link[index] = Node()
                current = current.link[index]
                
            self.insert_aux(current, key, i+1)

    def insert_2(self, key):
        """
        Insert the second key and its suffixes to the Trie
        by using a recursion. Used only after the first key is
        inserted.

        :Input:
        key: a string which will be inserted to the Trie along with
             its suffixes.

        :Output, return or postcondition: The key is inserted to the string
                                          along with its suffixes.

        :Time complexity: O(M^2), where M is the length of key
        :Aux space complexity: O(M^2)
        :Citation: FIT2004 2020sem02 Lecture11 Trie
        """
        maximum = 0

        for i in range(len(key) + 1):

            path = []

            current = self.root

            counter = 0

            counter = self.insert_aux_2(current, key, i, counter, path)

            # to determine the current longest non-branching
            # substring from the root. This will be stored in self.path
            if counter > maximum and len(path) > len(self.path):
                maximum = counter
                self.path = path

    def insert_aux_2(self, current, key, i, counter, path) -> None:
        """
        Recursion to insert the key[i::] into the Trie. Will return the
        path when inserting the current key.

        :Input:
        current : a node object which is the self.root
        key     : the string which is inserted into the Trie
        i       : an integer to determine which index of the key will be inserted
        counter : the length of the non-branching substring from the root (an integer)
        path    : an list of integers which determines the path of the substrings
                  from the root

        :Output, return or postcondition: The key[i::] is inserted to the string. Returns
                                          the how deep the substring from the root
                                          before it branches.
                                          

        :Time complexity: O(M), where M is the length of key
        :Aux space complexity: O(M)
        :Citation: FIT2004 2020sem02 Lecture11 Trie
        """
        if i == len(key):
            if current.link[0] is not None:
                current = current.link[0]
                current.branch = True

            else:
                current.link[0] = Node(2)
                current = current.link[0]
            
            return counter

        else: # i != len(key)
            index = self._get_index(key[i])

            # If the character already exist in the node
            if current.link[index] is not None:
                # If the substring is from the same string,
                # don't include it in the path
                if current.link[index].type != 2:
                    path.append(index)
                current = current.link[index]
                counter += 1

            else:
                current.link[index] = Node(2)
                current = current.link[index]

            counter = self.insert_aux_2(current, key, i+1, counter, path)

        return counter

def round(x):
    """
    Takes in a number and rounds in to the nearest integer

    :Input:
    x: a number which will be rounded to the nearest integer

    :Output, return or postcondition: Returns a rounded integer

    :Time complexity: O(1)
    :Aux space complexity: O(1)
    """
    y = x % 1
    # Rounded up
    if y >= 0.5:
        return int(x - y + 1)
    # Rounded down
    return int(x - y)

def compare_subs(submission1, submission2):
    """
    Compares the two submission and finds the longest matching
    substring. This function uses a suffix Trie. The submissions
    are inserted one by one into the Trie and it will
    return the similar substring and also the percentage of the substring
    in each submission.

    :Input:
    submission1: a string that may contain lowercase letters and spaces
    submission2: a string that may contain lowercase letters and spaces

    :Output, return or postcondition: Returns the longest matching substring
                                      between submission1 and submission2, and also
                                      the percentage of that substring in submission1
                                      and to respectively.

    :Time complexity: O(N^2 + M^2), where N is len(submission1) and M is len(submission2) and
                      is dominated by O(max(N^2, M^2))
    :Aux space complexity: O(max(N^2, M^2))
    """
    # O(1)
    similarity_detector = SuffixTrie()
    # O(N^2)
    similarity_detector.insert(submission1)
    # O(M^2)
    similarity_detector.insert_2(submission2)

    sentence = ''

    # O(min(N,M)) worst case where the shorter submission
    # matches all the substrings in the longer submission
    for index in similarity_detector.path:
        # index == 1 is space
        if index == 1:
            sentence += ' '
        else:
            sentence += chr(index+ 97 - 2)

    # If either submission has 0 length, then return empty string
    if len(submission1) == 0 or len(submission2) == 0:
        return ['', 0, 0]

    # count the respective percentage
    percentage_1 = round((len(sentence)/len(submission1)) * 100)
    percentage_2 = round((len(sentence)/len(submission2)) * 100)

    return [sentence, percentage_1, percentage_2]

