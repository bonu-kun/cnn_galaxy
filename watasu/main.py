from astropy.table import Table
import astropy.wcs
import astropy.io.fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from pydl.photoop.photoobj import unwrap_objid
import numpy as np
import matplotlib.pyplot as plt
import csv

DEG_CORRESPONDING_TO_5S = 0.00138889

def main():
    gz_table = Table.read('GalaxyZoo1_DR_table2.fits')

    # Galaxy Zooの表の上から 1000行目までを順に取り出す
    for galaxy_info in gz_table[:1000]:

        # Galaxy Zooの表からobjid, ra, decを取り出す
        objid = galaxy_info[0]
        ra_hms = galaxy_info[1]
        dec_hms = galaxy_info[2]

        # ra,decを60進数から10進数に変換
        c = SkyCoord(ra_hms + ' ' + dec_hms, unit=(u.hourangle, u.deg))
        ra_deg = c.ra.degree
        dec_deg = c.dec.degree

        # objidからfitsファイルを特定する情報を取り出す
        params = unwrap_objid(objid)
        run = str(params['run'])
        camcol = str(params['camcol'])
        frame = str(params['frame'])

        # objidの銀河が含まれるfitsファイルを開く
        fits = astropy.io.fits.open('../' + run.zfill(6) + '/fpC-' + run.zfill(6) + '-r' + camcol + '-' + frame.zfill(4) + '.fit.gz')
        data = fits[0].data
        header = fits[0].header
        
        # fitsファイルのheaderと銀河のra, decから,
        # 画像内の座標を得る
        wcs = astropy.wcs.WCS(header)
        px, py = wcs.wcs_world2pix(ra_deg, dec_deg, 0, ra_dec_order=True)
        px = int(px)
        py = int(py)
        

        # 得られた座標を中心として画像を切り出し、保存する
        r = 30
        img_galaxy = data[py-r:py+r, px-r:px+r]
        plt.title('RA = {}, Dec = {}'.format(ra_hms, dec_hms))
        plt.imshow(np.log10(img_galaxy.T[::-1, ::-1]), cmap='gray', vmin=3.05, vmax=3.09)

        # 画像上で5秒角に対応する長さの線を描画する
        if header['CTYPE1'] == 'DEC--TAN':
            pix_corresponding_to_5s = round(DEG_CORRESPONDING_TO_5S / header['cd1_1'])
        else:
            pix_corresponding_to_5s = round(DEG_CORRESPONDING_TO_5S / header['cd1_2'])
            
        # cd1_2 = abs(header['CD1_2'])
        # pix_corresponding_to_5s = round(DEG_CORRESPONDING_TO_5S / cd1_2)
        
        plt.hlines(y=17, xmin=10, xmax=10+pix_corresponding_to_5s)
        plt.savefig('data/img/' + str(objid), bbox_inches='tight')
        plt.close()

        
        # 画像をnpyで保存する
        np.save('data/npz/' + str(objid), img_galaxy)

        
        # 画像とGalaxy Zooのデータをnpzで保存する
        #np.savez('data/npz/' + str(objid), img_galaxy, galaxy_info)
        # img = Image.fromarray(np.uint8(fits[0].data[py-r:py+r, px-r:px+r]))
        
    csv_obj.close()
    
if __name__ == '__main__':
    main()
