import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_by_perlocal(reachDischargeFile, reachMetricsFile,gwDataFile, figFile, plotCol='Per_Local', axisTitle="Percent of baseflow discharged locally"):
    dischargeDF = pd.read_csv(reachDischargeFile)
    metricsDF = pd.read_csv(reachMetricsFile)
    gwData = np.load(gwDataFile)
    
    dischargeDF=dischargeDF.merge(metricsDF)
    dischargeDF = dischargeDF.loc[dischargeDF.variable=="temp"]
    
    gwTst = pd.DataFrame(gwData['GW_tst'], columns=gwData['GW_cols'])
    gwTst['partition']="tst"
    gwTrn = pd.DataFrame(gwData['GW_trn'], columns=gwData['GW_cols'])
    gwTrn['partition']="trn"
    gwVal = pd.DataFrame(gwData['GW_val'], columns=gwData['GW_cols'])
    gwVal['partition']="val"
    gwDF = gwTst.append(gwTrn,ignore_index=True).append(gwVal,ignore_index=True)
    
    gwDF['group']="Atmosphere"
    gwDF.loc[gwDF.Ar_obs.isnull(),'group']='Unknown'
    gwDF.loc[gwDF.delPhi_obs>=10,"group"]="Shallow"
    gwDF.loc[(gwDF.delPhi_obs<=10) & (gwDF.Ar_obs<0.65),"group"]="Deep"
    
    dischargeDF=dischargeDF.merge(gwDF)
    dischargeDF.replace([np.inf, -np.inf], np.nan, inplace=True)
    dischargeDF.loc[dischargeDF.nse.isnull(),"rmse"]=np.nan
    
    fig = plt.figure(figsize=(30, 15))
    thisFig = 0
    for thisPart in np.unique(dischargeDF.partition):
            thisData = dischargeDF.loc[dischargeDF.partition==thisPart]

            thisFig = thisFig + 1
            ax = fig.add_subplot(len(np.unique(dischargeDF.partition)), 4, thisFig)
            ax.set_title('RMSE versus Simulated Discharge, {}'.format(thisPart))
            #ax.axline((np.nanmean(thisData['{}_pred'.format(thisMetric)]),np.nanmean(thisData['{}_pred'.format(thisMetric)])), slope=1.0,linewidth=1, color='black', label="1 to 1 line")
            colorDict = {"Atmosphere":"red","Shallow":"green","Deep":"blue", "Unknown":"gray"}
            for thisGroup in np.unique(thisData['group'])[::-1]:
                thisColor = colorDict[thisGroup]
                ax.scatter(x=thisData.loc[thisData.group==thisGroup,plotCol],y=thisData.loc[thisData.group==thisGroup,"rmse"],label="RGCN - %s"%thisGroup,color=thisColor)

    #                ax.scatter(x=thisData['{}_obs'.format(thisMetric)],y=thisData['{}_sntemp'.format(thisMetric)],label="SNTEMP",color="red")
    #         for i, label in enumerate(thisData.seg_id_nat):
    #             ax.annotate(int(label), (thisData.loc[thisData.group==thisGroup,'Per_Local'][i],thisData.loc[thisData.group==thisGroup,"rmse"][i]))
            if thisFig==1:
                      ax.legend()
            ax.set_xlabel(axisTitle)
            ax.set_ylabel("RMSE, in degrees C")



            thisFig = thisFig + 1
            #colsToPlot = ['{}_obs'.format(thisMetric),'{}_pbm'.format(thisMetric),'{}_pred'.format(thisMetric)]
            nObs =["n: " + str(np.sum(np.isfinite(thisData.loc[thisData.group==x,plotCol].values))) for x in np.unique(thisData.group)]
            ax = fig.add_subplot(len(np.unique(dischargeDF.partition)), 4, thisFig)
            ax.set_title('{}, {}'.format(axisTitle,thisPart))
            ax=sns.boxplot(x='group',y=plotCol, data=thisData, order=np.unique(thisData.group), palette=colorDict)
            # Add it to the plot
            pos = range(len(nObs))
            for tick,label in zip(pos,ax.get_xticklabels()):
                ax.text(pos[tick],
                        np.nanmin(thisData[plotCol].values)-0.1*(np.nanmax(thisData[plotCol].values)-np.nanmin(thisData[plotCol].values)),
                        nObs[tick],
                        horizontalalignment='center',
                        weight='semibold')
            ax.set_ylim(np.nanmin(thisData[plotCol].values)-0.2*(np.nanmax(thisData[plotCol].values)-np.nanmin(thisData[plotCol].values)),np.nanmax(thisData[plotCol].values))

            thisFig = thisFig + 1
            ax = fig.add_subplot(len(np.unique(dischargeDF.partition)), 4, thisFig)
            ax.set_title('Ar versus Simulated Discharge, {}'.format(thisPart))
            #ax.axline((np.nanmean(thisData['{}_pred'.format(thisMetric)]),np.nanmean(thisData['{}_pred'.format(thisMetric)])), slope=1.0,linewidth=1, color='black', label="1 to 1 line")
            colorDict = {"Atmosphere":"red","Shallow":"green","Deep":"blue", "Unknown":"gray"}
            for thisGroup in np.unique(thisData['group'])[::-1]:
                thisColor = colorDict[thisGroup]
                ax.scatter(x=thisData.loc[thisData.group==thisGroup,plotCol],y=thisData.loc[thisData.group==thisGroup,"Ar_obs"],label="RGCN - %s"%thisGroup,color=thisColor)

    #                ax.scatter(x=thisData['{}_obs'.format(thisMetric)],y=thisData['{}_sntemp'.format(thisMetric)],label="SNTEMP",color="red")
    #         for i, label in enumerate(thisData.seg_id_nat):
    #             ax.annotate(int(label), (thisData.loc[thisData.group==thisGroup,'Per_Local'][i],thisData.loc[thisData.group==thisGroup,"rmse"][i]))
            if thisFig==1:
                      ax.legend()
            ax.set_xlabel(axisTitle)
            ax.set_ylabel("Ar, dimensionless")
            
            thisFig = thisFig + 1
            ax = fig.add_subplot(len(np.unique(dischargeDF.partition)), 4, thisFig)
            ax.set_title('delta Phi versus Simulated Discharge, {}'.format(thisPart))
            #ax.axline((np.nanmean(thisData['{}_pred'.format(thisMetric)]),np.nanmean(thisData['{}_pred'.format(thisMetric)])), slope=1.0,linewidth=1, color='black', label="1 to 1 line")
            colorDict = {"Atmosphere":"red","Shallow":"green","Deep":"blue", "Unknown":"gray"}
            for thisGroup in np.unique(thisData['group'])[::-1]:
                thisColor = colorDict[thisGroup]
                ax.scatter(x=thisData.loc[thisData.group==thisGroup,plotCol],y=thisData.loc[thisData.group==thisGroup,"delPhi_obs"],label="RGCN - %s"%thisGroup,color=thisColor)

    #                ax.scatter(x=thisData['{}_obs'.format(thisMetric)],y=thisData['{}_sntemp'.format(thisMetric)],label="SNTEMP",color="red")
    #         for i, label in enumerate(thisData.seg_id_nat):
    #             ax.annotate(int(label), (thisData.loc[thisData.group==thisGroup,'Per_Local'][i],thisData.loc[thisData.group==thisGroup,"rmse"][i]))
            if thisFig==1:
                      ax.legend()
            ax.set_xlabel(axisTitle)
            ax.set_ylabel("Ar, dimensionless")
            
            
    plt.savefig(figFile)