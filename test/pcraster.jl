using Conda
using PyCall
# Conda.add("pcraster"; channel="conda-forge")

pcr = pyimport("pcraster")
np = pyimport("numpy")

function spread_pcr(points, initial, friction; res=2, nodata=NaN)
    pcr.setclone(size(points)..., res, 0, 0)
    ppoints = pcr.numpy2pcr(pcr.Nominal, points, nodata)
    @time presult = pcr.spread(ppoints, initial, friction)
    transpose(pcr.pcr2numpy(presult, nodata))
end
