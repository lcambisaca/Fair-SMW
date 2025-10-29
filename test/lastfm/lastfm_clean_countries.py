import pandas as pd

df = pd.read_csv('lastfm_asia_target.csv')
# analyze the count of users in each country
print(df['target'].value_counts())
# delete users from countries that have less than 500 users
df.drop(df.loc[df['target'].isin([8, 5, 15, 16, 11, 7, 2, 13, 9, 12, 1, 4])].index, inplace=True)
print(df)
df.to_csv('lastfm_clean_countries.csv', header=False, index=False)
