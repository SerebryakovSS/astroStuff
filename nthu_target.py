#!/usr/bin/python3
# author: Serebryakov S., 2019
#https://astro.swarthmore.edu/transits.cgi
import os
import json
import requests
from lxml import html
page = requests.get('https://astro.swarthmore.edu/print_transits.cgi?observatory_string=28.758333%3B-17.88%3BAtlantic%2FCanary%3BRoque+de+los+Muchachos%2C+La+Palma&use_utc=1&observatory_latitude=&observatory_longitude=&timezone=UTC&start_date=today&days_to_print=1&days_in_past=0&minimum_start_elevation=+&and_vs_or=or&minimum_end_elevation=+&minimum_depth=0&target_string=&single_object=0&ra=&dec=&epoch=&period=&duration=&target=&show_ephemeris=0&print_html=1&twilight=-12&max_airmass=2.4')
document = html.fromstring(page.content)
data_arr  = []
date_par = ""
for tr_elements in document.xpath('//tbody/tr'):
    td_elements  = tr_elements.xpath('.//td/text() | .//td/span/text()' )
    target_names = tr_elements.xpath('.//td//span/a/text()')
    tt_magnitude = tr_elements.xpath('.//td//div/div/text()')[0]
    obs_end_time = tr_elements.xpath('.//td//span/span/text()')[0].replace('\n','')
    for i in range(len(td_elements)):
        td_elements[i] = td_elements[i].replace('\n','')
        td_elements[i] = td_elements[i].replace('\t','')
        td_elements[i] = td_elements[i].strip()
    date_par = td_elements.pop(0).split()[1]
    for target in target_names:
        data_arr.append({"date"           : date_par ,
                         "target"         : target   ,
                         "mag"            : tt_magnitude,
                         "start-end"      : str(td_elements[18])+" - "+str(obs_end_time),
                         "dur"            : td_elements[26],
                         "elev"           : [ elem.replace('°','') for elem in [td_elements[31],td_elements[34],td_elements[37]] ],
                         "RA,DEC"         : [ td_elements[-5],td_elements[-4] ],
                         "depth"          : td_elements[-1] })
print("=======================================================================================================================================================================")
target_list = []
for target in data_arr:
    if float(target["mag"]) < 16:
        dec = int(target["RA,DEC"][1][0:3])
        if dec > -30 and dec < 35:
            if int(target["elev"][0]) > 45 and int(target["elev"][1]) > 45:
                target_list.append(target)
if len(target_list) == 0:
    print("couldn't find any target for today...")
else:
    print("||%15s||%15s||%15s||%15s||%15s||%25s||%35s||%15s||"%("date","target","magnitude","start-end","duration","elev","RA,DEC","depth"))
    with open("trastits_tonight","w") as t_file:
        for target in target_list:
            t_file.write(json.dumps(target)+"\n")
            print("||%15s||%15s||%15s||%15s||%15s||%25s||%35s||%15s||"%(target['date'],target['target'],target['mag'],target['start-end'],target['dur'],target['elev'],target['RA,DEC'],target['depth']))
print("=======================================================================================================================================================================")    
'''
OUTPUT:
=======================================================================================================================================================================
||           date||         target||      magnitude||      start-end||       duration||                     elev||                             RA,DEC||          depth||
||     2019-08-05||      GJ 9827 c||          10.10||  02:34 - 04:24||           1:49||       ['56', '60', '58']||     ['23:27:04.62', '-01:17:12.5']||            0.4||
||     2019-08-05||       K2-141 b||          11.45||  03:13 - 04:09||           0:56||       ['60', '60', '59']||     ['23:23:39.97', '-01:11:21.4']||            0.4||
=======================================================================================================================================================================
'''
