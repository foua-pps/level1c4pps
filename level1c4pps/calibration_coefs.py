import datetime

CALIB_MODE = 'Nominal'
COEFS_MEIRINK = dict(
    MSG1=dict(
        VIS006=dict(b=24.346, a=0.0003739),
        VIS008=dict(b=30.989, a=0.0003111),
        IR_016=dict(b=22.869, a=0.0000065)
    ),
    MSG2=dict(
        VIS006=dict(b=21.026, a=0.0002556),
        VIS008=dict(b=26.875, a=0.0001835),
        IR_016=dict(b=21.394, a=0.0000498)
    ),
    MSG3=dict(
        VIS006=dict(b=19.829, a=0.0005856),
        VIS008=dict(b=25.284, a=0.0006787),
        IR_016=dict(b=23.066, a=-0.0000286)
    ),
    MSG4=dict(
        VIS006=dict(b=22.960, a=0),
        VIS008=dict(b=28.885, a=0),
        IR_016=dict(b=22.151, a=0)
    )
)

REF_DATE = datetime.date(2000, 1, 1)
REF_TIME = datetime.datetime(2000, 1, 1, 0, 0)

def calib_meirink(platform, channel, time):
    """Get MODIS-intercalibrated gain and offset for SEVIRI VIS channels

    Reference: http://msgcpp.knmi.nl/mediawiki/index.php/MSG-SEVIRI_solar_channel_calibration

    :returns: gain, offset [mW m-2 sr-1 (cm-1)-1]
    """
    if time < REF_TIME:
        raise ValueError('Given time ({0}) is < reference time ({1})'.format(time, REF_TIME))

    a = COEFS_MEIRINK[platform][channel]['a']
    b = COEFS_MEIRINK[platform][channel]['b']
    delta_days = (time - REF_TIME).total_seconds() / 3600.0 / 24.0
    gain = (b + a * delta_days) / 1000.0  # micro Watts -> milli Watts
    offset = -51.0 * gain  # Space count is 51

    return gain, offset

def calib_meirink_date(platform, channel, date):
    """Get MODIS-intercalibrated gain and offset for SEVIRI VIS channels

    Reference: http://msgcpp.knmi.nl/mediawiki/index.php/MSG-SEVIRI_solar_channel_calibration

    :returns: gain, offset [mW m-2 sr-1 (cm-1)-1]
    """
    if date < REF_DATE:
        raise ValueError('Given date ({0}) is < reference date ({1})'.format(date, REF_DATE))

    a = COEFS_MEIRINK[platform][channel]['a']
    b = COEFS_MEIRINK[platform][channel]['b']
    gain = (b + a*(date - REF_DATE).days) / 1000.0  # micro Watts -> milli Watts
    offset = -51.0 * gain  # Space count is 51

    return gain, offset

def get_calibration_for_time(platform, time):
    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink(platform=platform, channel=channel, time=time)
        coefs[channel] =  {'gain': gain, 'offset': offset}

    return coefs 

def get_calibration_for_date(platform, date):
    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink_date(platform=platform, channel=channel, date=date)
        coefs[channel] =  {'gain': gain, 'offset': offset}

    return coefs

def test_get_calibration_for_date():
    coefs = get_calibration_for_date(platform='MSG3', 
                                     date=datetime.date(2018, 1, 18))
    REF = {'VIS006': {'gain': 0.023689275200000002, 'offset': -1.2081530352}, 
           'VIS008': {'gain': 0.029757990399999996, 'offset': -1.5176575103999999},
           'IR_016': {'gain': 0.0228774688, 'offset': -1.1667509087999999}}
    for channel in REF.keys():
        if (REF[channel]['gain'] == coefs[channel]['gain'] and
            REF[channel]['offset'] == coefs[channel]['offset']):
            print("Calibration for channel {:s} is OK!".format(channel))
                  

def test_get_calibration_for_time():
    coefs = get_calibration_for_time(platform='MSG3', 
                                     time=datetime.datetime(2018, 1, 18, 0, 0))
    REF = {'VIS006': {'gain': 0.023689275200000002, 'offset': -1.2081530352}, 
           'VIS008': {'gain': 0.029757990399999996, 'offset': -1.5176575103999999},
           'IR_016': {'gain': 0.0228774688, 'offset': -1.1667509087999999}}
    for channel in REF.keys():
        if (REF[channel]['gain'] == coefs[channel]['gain'] and
            REF[channel]['offset'] == coefs[channel]['offset']):
            print("Calibration for channel {:s} is OK!".format(channel))
        else :
            print(REF[channel]['gain'] - coefs[channel]['gain'])
        
def test_get_calibration():
    coefs1 = get_calibration_for_time(platform='MSG3', 
                                      time=datetime.datetime(2018, 1, 18, 23, 59))
    coefs2 = get_calibration_for_date(platform='MSG3', 
                                      date=datetime.date(2018, 1, 19))
    for channel in coefs1.keys():
        if (coefs1[channel]['gain'] - coefs2[channel]['gain'] < 10e-8 and
            coefs1[channel]['offset'] - coefs2[channel]['offset'] < 10e-8):
            print("Calibration for channel {:s} is OK!".format(channel))
        else :
            print(coefs2[channel]['gain'] - coefs1[channel]['gain'])
                  

if __name__ == '__main__':
    test_get_calibration_for_time()
    test_get_calibration_for_date()
    test_get_calibration()
    time = datetime.datetime(2018, 1, 18, 12, 0)
    platform = 'MSG3'

    coefs = {}
    for channel in ('VIS006', 'VIS008', 'IR_016'):
        gain, offset = calib_meirink(platform=platform, channel=channel, time=time)
        coefs[channel] = {'gain': gain, 'offset': offset}

    print(coefs)


