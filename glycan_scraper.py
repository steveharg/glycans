import urllib2
from bs4 import BeautifulSoup
import pickle

sample_type = 'cell'
# sample_type = 'tissue'

glypage_home = 'http://www.functionalglycomics.org'

if sample_type == 'cell':
    msa_files_dir = 'cell_type_msa_files'
    gly_page = 'http://www.functionalglycomics.org/glycomics/common/jsp/samples/searchSample.jsp?templateKey=2&12=CellType&operation=refine'
else:
    msa_files_dir = 'tissue_type_msa_files'
    gly_page = 'http://www.functionalglycomics.org/glycomics/common/jsp/samples/searchSample.jsp?templateKey=1&12=Tissue&operation=refine'



page = urllib2.urlopen(gly_page)
soup = BeautifulSoup(page, 'html.parser')

gtTablePresentation_trs = soup.find_all('tr', attrs={'class': 'gtTablePresentation'})

i = 0

for gtTablePresentation_tr in gtTablePresentation_trs:
    glycoEnzymeTablePresentation_tds = gtTablePresentation_tr.find_all('td', attrs={'class': 'glycoEnzymeTablePresentation'})
    for glycoEnzymeTablePresentation_td in glycoEnzymeTablePresentation_tds:
        anchors = glycoEnzymeTablePresentation_td.find_all('a')
        for anchor in anchors:
            imgs = anchor.find_all('img')
            for img in imgs:
                img_attrs = img.attrs
                if 'src' in img_attrs:
                    if 'MSA.jpg' in img_attrs['src']:
                        print(img)
                        msa_url = anchor.get('href')

                        webSiteBodyText_tds = gtTablePresentation_tr.find_all('td', attrs={'class': 'webSiteBodyText'})

                        sample_metadata = {}

                        if sample_type == 'cell':
                            metadata_keys = ['species', 'comments', 'cell_type', 'glycan_type_analyzed',
                                             'participating_investigator']
                        else:
                            metadata_keys = ['mice_clony_code', 'species', 'mouse_strain', 'tissue',
                                             'comments', 'glycan_type_analyzed', 'mouse_type', 'patient_code',
                                             'participating_investigator']

                        metadata_idx = 0
                        for metadata in webSiteBodyText_tds:
                            sample_metadata[metadata_keys[metadata_idx]] = metadata.text.strip()
                            metadata_idx += 1

                        msa_page = urllib2.urlopen(glypage_home+'/'+msa_url)
                        msa_page_soup = BeautifulSoup(msa_page, 'html.parser')

                        msa_base_file_path =  '/coreCStatic/allmsdfiles/'
                        msa_url_split = msa_url.split(msa_base_file_path)
                        msa_file_name = msa_url_split[1]

                        with open(msa_files_dir + '/' + msa_file_name, 'w') as msa_file:
                            msa_file.write(msa_page_soup.text)

                        pickle.dump(sample_metadata, open(msa_files_dir + "/" + msa_file_name + ".pickle", "wb"))

    i += 1

print(i)



