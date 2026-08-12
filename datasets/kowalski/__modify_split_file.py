import pandas as pd

split_df = pd.read_csv('split_trainonhek293ft_r2.csv')
print(split_df.groupby(['split', 'subsplit']).size())
split_df = split_df.query('subsplit != "HEK293FT_R2" or split != "test"')
print(split_df.groupby(['split', 'subsplit']).size())
split_df['subsplit'] = split_df['subsplit'].str.split('_').str[0]
print(split_df.groupby(['split', 'subsplit']).size())
split_df.to_csv('split_trainonhek293ft_r2_modified.csv', index=False)


split_df = pd.read_csv('split_trainonk562_r1.csv')
print(split_df.groupby(['split', 'subsplit']).size())
split_df = split_df.query('subsplit != "K562_R1" or split != "test"')
print(split_df.groupby(['split', 'subsplit']).size())
split_df['subsplit'] = split_df['subsplit'].str.split('_').str[0]
print(split_df.groupby(['split', 'subsplit']).size())
split_df.to_csv('split_trainonk562_r1_modified.csv', index=False)

