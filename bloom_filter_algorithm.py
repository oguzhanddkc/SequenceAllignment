# Python 3 program to build Bloom Filter
# Install mmh3 and bitarray 3rd party module first
# pip install mmh3
# pip install bitarray
import math
import mmh3
from bitarray import bitarray
from random import shuffle
import kmer_algorithm
import sqlite3

class BloomFilter(object):
    '''
    Class for Bloom filter, using murmur3 hash function
    '''

    def __init__(self, items_count, fp_prob):
        '''
        items_count : int
            Number of items expected to be stored in bloom filter
        fp_prob : float
            False Positive probability in decimal
        '''
        # False possible probability in decimal
        self.fp_prob = fp_prob

        # Size of bit array to use
        self.size = self.get_size(items_count, fp_prob)

        # number of hash functions to use
        self.hash_count = self.get_hash_count(self.size, items_count)

        # Bit array of given size
        self.bit_array = bitarray(self.size)

        # initialize all bits as 0
        self.bit_array.setall(0)

    def add(self, item):
        '''
        Add an item in the filter
        '''
        digests = []
        for i in range(self.hash_count):
            # create digest for given item.
            # i work as seed to mmh3.hash() function
            # With different seed, digest created is different
            digest = mmh3.hash(item, i) % self.size
            digests.append(digest)

            # set the bit True in bit_array
            self.bit_array[digest] = True

    def check(self, item):
        '''
        Check for existence of an item in filter
        '''
        for i in range(self.hash_count):
            digest = mmh3.hash(item, i) % self.size
            if self.bit_array[digest] == False:
                # if any of bit is False then,its not present
                # in filter
                # else there is probability that it exist
                return False
        return True

    @classmethod
    def get_size(self, n, p):
        '''
        Return the size of bit array(m) to used using
        following formula
        m = -(n * lg(p)) / (lg(2)^2)
        n : int
            number of items expected to be stored in filter
        p : float
            False Positive probability in decimal
        '''
        m = -(n * math.log(p)) / (math.log(2) ** 2)
        return int(m)

    @classmethod
    def get_hash_count(self, m, n):
        '''
        Return the hash function(k) to be used using
        following formula
        k = (m/n) * lg(2)

        m : int
            size of bit array
        n : int
            number of items expected to be stored in filter
        '''
        k = (m / n) * math.log(2)
        return int(k)

    @staticmethod
    def run_bloom_filter(text, patterns):
        n = 20  # no of items to add
        p = 0.05  # false positive probability

        bloomf = BloomFilter(n, p)

        """print("Size of bit array:{}".format(bloomf.size))
        print("False positive Probability:{}".format(bloomf.fp_prob))
        print("Number of hash functions:{}".format(bloomf.hash_count))"""

        #text = "CTCAAAAGTCTGTGCCCCACCCACAGCGCCAGTGCACACAGCTTTCT"

        text_size, parca_size = len(text), 10
        text_parcalar = [text[i:i + parca_size] for i in range(0, text_size, parca_size)]

        """
        #word_present bizim textimizin 10luk parçalara ayırılmış halidir
        word_present = ["CTCAAAAGTC",
        				"TGTGCCCCAC",
        				"CCACAGCGCC",
        				"TGCACACAGCT"]
        """
        """
        #bu bizim patternler listemiz
        word_absent = ["CCACCCCGCGGGTCGAACA",
                    "CTCAAAAGTC",
                    "CGAGTTTGGGGCATCATA",
                    "GAGCCCAAGGCGATCTTTCC",
                    "GCCCGCAGGGCTTCCTGGT",
                    "GGACACTTTTGGGAGCGTG",
                    "CTCAAAAGTCT",
                    "aleyna"]
        """
        # parçalanmış textleri bloom filtera atar
        for item in text_parcalar:
            bloomf.add(item)

        shuffle(text_parcalar)
        shuffle(patterns)

        test_words = text_parcalar[:10] + patterns
        shuffle(test_words)
        #print(test_words)
        conn = sqlite3.connect("database.db")
        cursor = conn.cursor()
        sqlcommand = """
               CREATE TABLE IF NOT EXISTS BLOOMFILTER
               (
               pattern TEXT PRIMARY KEY,
               output TEXT
               )      
               """
        cursor = conn.execute(sqlcommand)

        """
        for word in test_words:
            if bloomf.check(word):
                if word in patterns:
                    print("'{}' is a false positive!".format(word))
                else:
                    print("'{}' is probably present!".format(word))
            else:
                print("'{}' is definitely not present!".format(word))
        """

        return test_words

    def insertVaribleIntoTable(text, start):
        try:
            sqliteConnection = sqlite3.connect('database.db')
            cursor = sqliteConnection.cursor()
            print("Connected to SQLite")

            sqlite_insert_with_param = """INSERT INTO BLOOMFILTER
                              (pattern, output) 
                              VALUES (?, ?);"""

            data_tuple = (text, start)
            cursor.execute(sqlite_insert_with_param, data_tuple)
            sqliteConnection.commit()
            print("Python Variables inserted successfully into SqliteDb_developers table")

            cursor.close()

        except sqlite3.Error as error:
            print("Failed to insert Python variable into sqlite table", error)
        finally:
            if sqliteConnection:
                sqliteConnection.close()
                print("The SQLite connection is closed")


if __name__ == "__main__":
    n = 20  # no of items to add
    p = 0.05  # false positive probability
    bloomf = BloomFilter(n, p)
    bloomf.run_bloom_filter()
    """
    n = 20  # no of items to add
    p = 0.05  # false positive probability

    bloomf = BloomFilter(n, p)
    print("Size of bit array:{}".format(bloomf.size))
    print("False positive Probability:{}".format(bloomf.fp_prob))
    print("Number of hash functions:{}".format(bloomf.hash_count))

    # words to be added
    word_present = ['abound', 'abounds', 'abundance', 'abundant', 'accessible',
                    'bloom', 'blossom', 'bolster', 'bonny', 'bonus', 'bonuses',
                    'coherent', 'cohesive', 'colorful', 'comely', 'comfort',
                    'gems', 'generosity', 'generous', 'generously', 'genial']

    # word not added
    word_absent = ['bluff', 'cheater', 'hate', 'war', 'humanity',
                   'racism', 'hurt', 'nuke', 'gloomy', 'facebook',
                   'geeksforgeeks', 'twitter']

    for item in word_present:
        bloomf.add(item)

    shuffle(word_present)
    shuffle(word_absent)

    test_words = word_present[:10] + word_absent
    shuffle(test_words)
    for word in test_words:
        if bloomf.check(word):
            if word in word_absent:
                print("'{}' is a false positive!".format(word))
            else:
                print("'{}' is probably present!".format(word))
        else:
            print("'{}' is definitely not present!".format(word))
    """
