#!/usr/bin/python

def mysql_conn():

    import os, MySQLdb
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="ncbi_taxonomy") # name of the data base

    cursor = conn.cursor()

    return conn, cursor

def get_ncbi_taxon_data(taxon_id):

    conn, cursor = mysql_conn()

    sql = 'select ncbi_taxon_id,node_rank, name, parent_taxon_id, t1.taxon_id from taxon_name t1 ' \
          ' inner join taxon t2 on t1.taxon_id=t2.taxon_id' \
          ' where ncbi_taxon_id=%s and name_class="scientific name";' % taxon_id
    cursor.execute(sql)
    print sql
    data = cursor.fetchall()[0]
    return data

def get_sql_taxon_data(taxon_id):

    conn, cursor = mysql_conn()

    sql = 'select ncbi_taxon_id,node_rank, name, parent_taxon_id, t1.taxon_id from taxon_name t1 ' \
          ' inner join taxon t2 on t1.taxon_id=t2.taxon_id' \
          ' where t1.taxon_id=%s and name_class="scientific name";' % taxon_id
    cursor.execute(sql)
    data = cursor.fetchall()[0]
    return data


def get_child(taxon_id):

    conn, cursor = mysql_conn()

    sql = 'select ncbi_taxon_id,node_rank, name, parent_taxon_id, t1.taxon_id from taxon_name t1 ' \
          ' inner join taxon t2 on t1.taxon_id=t2.taxon_id where parent_taxon_id=%s;' % taxon_id
    cursor.execute(sql)
    data = cursor.fetchall()[0]
    return data

def get_taxonomy_downstream_path(taxon_id):
    ncbi_taxon_id , node_rank , name , parent_taxon_id, sql_taxon_id = get_ncbi_taxon_data(taxon_id)
    rank2data = {}
    rank2data[node_rank] = [ncbi_taxon_id , name , parent_taxon_id, sql_taxon_id]
    while sql_taxon_id != parent_taxon_id:
        ncbi_taxon_id, node_rank, name, parent_taxon_id, sql_taxon_id = get_sql_taxon_data(parent_taxon_id)
        rank2data[node_rank] = [ncbi_taxon_id, name, parent_taxon_id, sql_taxon_id]
    return rank2data

def get_taxonomy_upstream_path(taxon_id):
    ncbi_taxon_id , node_rank , name , parent_taxon_id, sql_taxon_id = get_ncbi_taxon_data(taxon_id)
    rank2data = {}
    rank2data[node_rank] = [ncbi_taxon_id , name , parent_taxon_id, sql_taxon_id]
    while sql_taxon_id != parent_taxon_id:
        ncbi_taxon_id, node_rank, name, parent_taxon_id, sql_taxon_id = get_sql_taxon_data(parent_taxon_id)
        rank2data[node_rank] = [ncbi_taxon_id, name, parent_taxon_id, sql_taxon_id]
    return rank2data

#print get_taxonomy_upstream_path("426430")['species']

