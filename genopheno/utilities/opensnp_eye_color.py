import pandas as pd


def eye_color_normalize(eye_color):
    # hazel is brown
    # things like green-brown count as green
    # things like brown-green
    if eye_color == '-':
        # Unknown. There is no reason to try to resolve this to blue_green or brown.
        return eye_color

    # find first index of blue_green or brown
    eye_color = eye_color.lower()
    blue_index = -1
    for blue_find in [eye_color.find('blue'), eye_color.find('green'), eye_color.find('hazel')]:
        if blue_find > -1:
            if blue_index == -1 or blue_find < blue_index:
                blue_index = blue_find

    brown_index = eye_color.find('brown')

    if brown_index > -1 and blue_index == -1:
        return 'Brown'
    elif blue_index > -1 and brown_index == -1:
        return 'Blue_Green'
    elif brown_index > -1 and blue_index > -1:
        # if the descriptions contains both blue_green and brown then choose the color that comes first
        if brown_index < blue_index:
            return 'Brown'
        else:
            return 'Blue_Green'
    else:
        # '-' represents unknown in OpenSNP
        return '-'


def eye_color_pheno():
    df = pd.read_csv('/Users/rob/Downloads/opensnp_datadump.current/phenotypes_201705311214.csv', error_bad_lines=False,
                     sep=';')

    dfnew = pd.DataFrame()
    dfnew['user_id'] = df['user_id']
    pheno_col = 'phenotype'
    dfnew[pheno_col] = df['Eye color']

    # normalize eye colors
    dfnew[pheno_col] = dfnew[pheno_col].apply(eye_color_normalize)

    # remove no responses and non standard responses
    dfnew.drop(dfnew[dfnew[pheno_col] == '-'].index, inplace=True)

    # save
    dfnew.to_csv('/Users/rob/Downloads/opensnp_datadump.current/eye_color.csv', index=False)


if __name__ == '__main__':
    eye_color_pheno()
