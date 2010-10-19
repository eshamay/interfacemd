import Utility
import ColumnPairPlotter as CPP

files = [i for i in Utility.SearchDirectoryTree('.','so2-bond*dat')]
CPP.PlotOffsetColumnPairs (files, 0)

