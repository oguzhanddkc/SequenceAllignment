import timeit
import sqlite3

EOS = "$"

txt_sequence = ""
ptrn_sequence = ""


def print_table(table):
    for i in table:
        print(i)


def print_dictionary(dictionary):
    for i in sorted(dictionary.keys()):
        print(str(i) + '\t' + str(dictionary[i]))


def bwt_table_construction(x):
    x = x + EOS
    bwtArray = []
    bwtArray.append(x)
    for i in range(1, len(x), 1):
        bwtArray.append(x[i::] + x[:i:])

    bwtArray.sort()
    return bwtArray


def bwt_string(array):
    output_string = ""
    for i in array:
        output_string += i[-1]
    return output_string


def sa_array(x):
    saDic = {}
    saArray = []
    for i in range(len(x) + 1):
        saDic[x[i::]] = i
    # Stupid merge sort
    for seq in sorted(saDic.keys()):
        saArray.append(saDic[seq])
    return saArray


def bwt_table_construction_equation(x, suffixArray):
    bwtString = []
    for i in range(len(x) + 1):
        if suffixArray[i] == 0:
            bwtString.append("$")
        else:
            bwtString.append(x[suffixArray[i] - 1])
    # bwtStringbwtString[i]
    return ''.join(bwtString)


def first_col(saArray, txt_sequence):
    x = ""
    txt_sequence = txt_sequence + EOS
    # print("sa array: {}".format(saArray))
    for i in saArray:
        x = x + str(txt_sequence[i])
    return x


saArray = sa_array(txt_sequence)
# print("Suffix Array: {}".format(str(saArray)))

bwtArray = bwt_table_construction(txt_sequence)
# print("raw bwt is: {}".format(bwtArray))

# print("FM Table:")
# print_table(bwtArray)

bwt_table_construction_equation(txt_sequence, saArray)
bwtString = bwt_string(bwtArray)


# print("**BWT string is {}".format(bwtString))


def count_table(bwtString):
    # counttalbe bir dict
    countTable = {}
    alphbet = sorted(set(bwtString))
    alphbetRev = alphbet
    alphbetRev.reverse()
    # ("bwt string'i oluşturan karakterlerin sıralanmış ve ters çevrilmiş hali: {}".format(alphbetRev))
    for i in alphbet:
        countTable[i] = 0
    # counts = 0
    # BWTString'de olan elemanların frequency'sini hesaplar
    for letter in bwtString:
        countTable[letter] += 1
    # print("count table: {}".format(countTable))
    bwtStringLetterCounts = len(bwtString)

    for letter in alphbetRev:
        # print("before: " + str(countTable[letter]))
        countTable[letter] = bwtStringLetterCounts - countTable[letter]
        # print("after: " + str(countTable[letter]))
        bwtStringLetterCounts = countTable[letter]

    return countTable


countTable = count_table(bwtString)


# print("count table:")
# print_dictionary(countTable)

def occ_table(bwtString):
    # dict
    occTable = {}
    alphbet = sorted(set(bwtString))
    # print alphbet
    for nt in alphbet:
        occTable[nt] = [0] * len(bwtString)
    # print("occtable: {}".format(occTable))
    for step in range(len(bwtString)):
        nt = bwtString[step]
        occTable[nt][step] = 1
    # print("occtable before: {}".format(occTable))
    for nt in alphbet:
        for step in range(1, len(bwtString), 1):
            occTable[nt][step] += occTable[nt][step - 1]
        # print("occtable after : {}".format(occTable))
    return occTable


occTable = occ_table(bwtString)
# print("occtable:")
# print_dictionary(occTable)


firstCol = first_col(saArray, txt_sequence)


# print("first column: {}".format(firstCol))

def letter_range_in_string(letter, string):
    # string = string[begin:end+1]
    # string = list(string)
    letterPositions = []
    letterPositions = [i for i, x in enumerate(list(string)) if x == letter]
    # return [min(letterPositions),max(letterPositions)]
    return letterPositions


def letter_LF(letter, oldposition, bwtString, saArray, countTable, occTable):
    # print firstCol[oldposition]
    # print bwtString[oldposition]
    # print countTable[bwtString[oldposition]]
    # print occTable[bwtString[oldposition]][oldposition]
    # newposition= countTable[bwtString[oldposition]]+occTable[bwtString[oldposition]][oldposition]-1
    newposition = countTable[letter] + occTable[letter][oldposition] - 1
    return newposition


def string_search(ptrn_sequence, firstCol, bwtString, saArray, countTable, occTable):
    hitSApositions = []
    ptrn_length = len(ptrn_sequence)
    ptrn_sequence = list(ptrn_sequence)
    initLetter = ptrn_sequence.pop()
    startOld = min(letter_range_in_string(initLetter, firstCol))
    endOld = max(letter_range_in_string(initLetter, firstCol))
    if ptrn_length == 1:
        hitSApositions = letter_range_in_string(initLetter, firstCol)
    else:
        i = 2
        while i <= ptrn_length:
            letter = ptrn_sequence.pop()
            startOld = letter_LF(letter, startOld, bwtString, saArray, countTable, occTable)
            endOld = letter_LF(letter, endOld, bwtString, saArray, countTable, occTable)
            i += 1
        hitSApositions = [startOld, endOld]
        for i in range(len(hitSApositions)):
            hitSApositions[i] = saArray[hitSApositions[i]]

    return list(set(hitSApositions))


def align(txt_sequence, ptrn_sequence):
    start_prog = run_time_start()

    saArray = sa_array(txt_sequence)
    bwtArray = bwt_table_construction(txt_sequence)
    bwt_table_construction_equation(txt_sequence, saArray)
    bwtString = bwt_string(bwtArray)
    countTable = count_table(bwtString)
    occTable = occ_table(bwtString)
    firstCol = first_col(saArray, txt_sequence)
    result = string_search(ptrn_sequence, firstCol, bwtString, saArray, countTable, occTable)
    # print("Your exact alignment is located at: {}".format(result))

    stop_prog = run_time_stop()
    run_time = stop_prog - start_prog

    return result, run_time


"""
def run_time():
    
    align(txt_sequence, ptrn_sequence)
    
    return run_time
"""


def run_time_start():
    start = timeit.default_timer()
    return start


def run_time_stop():
    stop = timeit.default_timer()
    return stop


conn = sqlite3.connect("database.db")
cursor = conn.cursor()
sqlcommand = """
    CREATE TABLE IF NOT EXISTS FMINDEX
    (
    pattern TEXT PRIMARY KEY,
    found_index TEXT
    )      
    """
cursor = conn.execute(sqlcommand)


def insertVaribleIntoTable(text, start):
    try:
        sqliteConnection = sqlite3.connect('database.db')
        cursor = sqliteConnection.cursor()
        print("Connected to SQLite")

        sqlite_insert_with_param = """INSERT INTO FMINDEX
                          (pattern, found_index) 
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
# son = align(txt_sequence, ptrn_sequence)
# print(son)
