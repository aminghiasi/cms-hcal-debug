import ROOT as r
import sys

r.gROOT.SetBatch()
r.gStyle.SetOptStat(0)

infile, outfile = sys.argv[1:]

f = r.TFile(infile)
tps = f.Get("compareL1T/tps")

c = r.TCanvas()
c.SetRightMargin(c.GetRightMargin() * 1.5)
c.SaveAs(outfile + '[')
tps.Draw("iphi:ieta>>histdet(85,-42.5,42.5,73,0,72)", "(version==1) * (fg1-fg1_emul)", "COLZ")
r.gDirectory.Get("histdet").SetTitle(";ieta;iphi;#sum L1T FG - FG reemulated")
c.SaveAs(outfile)
c.SetLogz()
tps.Draw("iphi:ieta>>histdet2(85,-42.5,42.5,73,0,72)", "(version==1) * ((fg1 && soi>0) && !fg1_emul)", "COLZ")
tps.Draw("iphi:ieta>>histdet3(85,-42.5,42.5,73,0,72)", "(version==1)", "COLZ")
r.gDirectory.Get("histdet2").Divide(r.gDirectory.Get("histdet3"))
r.gDirectory.Get("histdet2").SetTitle("FG set in RAW w/ SOI>0 but not re-emulated (rel. freq.);ieta;iphi;(#sum L1T FG && !FG reemulated) / (#sum # TP)")
r.gDirectory.Get("histdet2").GetZaxis().SetRangeUser(0.0001, 1)
r.gDirectory.Get("histdet2").Draw("COLZ")
c.SaveAs(outfile)
tps.Draw("iphi:ieta>>histdet2(85,-42.5,42.5,73,0,72)", "(version==1) * (!fg1 && fg1_emul)", "COLZ")
tps.Draw("iphi:ieta>>histdet3(85,-42.5,42.5,73,0,72)", "(version==1)", "COLZ")
r.gDirectory.Get("histdet2").Divide(r.gDirectory.Get("histdet3"))
r.gDirectory.Get("histdet2").SetTitle("FG re-emulated but not set in RAW (rel. freq.);ieta;iphi;(#sum !L1T FG && FG reemulated) / (#sum # TP)")
r.gDirectory.Get("histdet2").GetZaxis().SetRangeUser(0.0001, 1)
r.gDirectory.Get("histdet2").Draw("COLZ")
c.SaveAs(outfile)
tps.Draw("iphi:ieta>>histdet2(85,-42.5,42.5,73,0,72)", "(version==1) * fg1_emul", "COLZ")
tps.Draw("iphi:ieta>>histdet3(85,-42.5,42.5,73,0,72)", "(version==1)", "COLZ")
r.gDirectory.Get("histdet2").Divide(r.gDirectory.Get("histdet3"))
r.gDirectory.Get("histdet2").SetTitle("FG set when re-emulated (rel. freq.);ieta;iphi;(#sum FG reemulated) / (#sum # TP)")
r.gDirectory.Get("histdet2").GetZaxis().SetRangeUser(0.0001, 1)
r.gDirectory.Get("histdet2").Draw("COLZ")
c.SaveAs(outfile)
tps.Draw("fg1-fg1_emul:soi>>histsoi", "version==1", "COLZ")
r.gDirectory.Get("histsoi").SetTitle(";SOI;L1T FG - FG reemulated;")
c.SetLogz()
c.SaveAs(outfile)
c.SaveAs(outfile + ']')
