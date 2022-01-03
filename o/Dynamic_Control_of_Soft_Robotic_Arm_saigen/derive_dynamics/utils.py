"""いろいろ"""

import pickle

def save_eq_by_txt(eq, dir_path, file_name,):
    """txtで保存"""
    
    f = open(dir_path + file_name + ".txt", 'w')
    f.write(str(eq))
    f.close()


def save_obj_by_picke(obj, dir_path, file_name,):
    """バイナリで保存"""
    
    f = open(dir_path + file_name + ".binaryfile", 'wb')
    pickle.dump(obj, f)
    f.close()


def load_obj_from_picle(dir_path, file_name,):
    """バイナリからオブジェクトを復元"""
    
    f = open(dir_path + file_name + ".binaryfile", 'rb')
    return pickle.load(f)