import pandas as pd

df1 = pd.read_csv('lastfm_clean_countries.csv', header=None)
edge_list = df1[df1.columns[0]].tolist()
n = len(edge_list)
print(n)

df2 = pd.read_csv('lastfm_asia_edges.csv')
print(df2)

# keep edges that both nodes have more than 500 users
df3 = df2.loc[(df2['node_1'].isin(edge_list)) & (df2['node_2'].isin(edge_list))]
df3 = df3.reset_index()  # make sure indexes pair with number of rows
print(df3)

df = pd.DataFrame(columns=["node_1", "node_2"])
# change node indexing to 1 to n
for index, row in df3.iterrows():
    i = edge_list.index(row['node_1'])
    j = edge_list.index(row['node_2'])
    # print(i, j)
    entry = pd.DataFrame.from_dict({
        "node_1": [i],
        "node_2":  [j]
    })
    df = pd.concat([df, entry], ignore_index=True)

print(df)



# delete users from countries that have less than 500 users
# df.drop(df.loc[df['target'].isin([8, 5, 15, 16, 11, 7, 2, 13, 9, 12, 1, 4])].index, inplace=True)
# print(df)

df.to_csv('lastfm_clean_edges.csv', header=False, index=False)
