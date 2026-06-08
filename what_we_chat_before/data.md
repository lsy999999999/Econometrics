可以。我们这一步原始数据链接建议放进一个 `data_sources.txt`，因为老师要求压缩包里包括原始数据，数据集较大时可以用 `.txt` 给出数据链接地址。

## 1. 原始数据下载网址

**A. 49行业组合日度收益率**

Kenneth French Data Library 官方页面里有 “49 Industry Portfolios [Daily]” 的 TXT、CSV 和 Details 链接。([达特茅斯大学塔克商学院][1])

```txt
官方页面：
https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

直接CSV压缩包：
https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/49_Industry_Portfolios_daily_CSV.zip

49行业说明页面：
https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/det_49_ind_port.html
```

**B. 2年期美债收益率 DGS2**

FRED 的 DGS2 是 “Market Yield on U.S. Treasury Securities at 2-Year Constant Maturity”，单位是 percent，非季调，频率是 daily。([FRED][2])

```txt
官方页面：
https://fred.stlouisfed.org/series/DGS2

直接CSV：
https://fred.stlouisfed.org/graph/fredgraph.csv?id=DGS2
```

**C. FOMC事件日期与声明**

美联储 FOMC 日历页面提供会议日期、声明和会议纪要链接；美联储也说明 FOMC 通常每年有8次定期会议，必要时也会召开其他会议。([联邦储备委员会][3])

```txt
FOMC会议日历、声明、纪要：
https://www.federalreserve.gov/monetarypolicy/fomccalendars.htm

2022-11-02 FOMC声明：
https://www.federalreserve.gov/newsevents/pressreleases/monetary20221102a.htm

2022-11-02 Implementation Note：
https://www.federalreserve.gov/newsevents/pressreleases/monetary20221102a1.htm
```

## 2. 为什么先选 2022-11-02？

我选 **2022-11-02** 不是说它是唯一正确日期，而是把它作为“第一步简单截面回归”的示范事件。这个日期有几个优点。

第一，它是一个方向非常明确的紧缩性 FOMC 事件。美联储当天声明明确说，将联邦基金利率目标区间提高到 **3.75%–4.00%**，并且认为继续加息是适当的，以使货币政策立场足够限制性，从而让通胀回到2%。([联邦储备委员会][4]) 这和我们的假设 `H1: β1 < 0` 是匹配的：如果事件是紧缩性的，那么更怕利率上升的行业应该受到更大负面冲击。

第二，这个事件处在 2022 年快速加息周期中，市场对利率、贴现率和成长股估值都非常敏感，因此用行业对 DGS2 的历史敏感度解释事件窗口收益，有比较清楚的资产定价含义。DGS2 本身是日度、百分比、非季调数据，适合和日度行业收益匹配。([FRED][2])

第三，它是正常 FOMC 会议和正常声明，不是疫情期间的紧急会议，也不是零利率下限附近的非常规政策事件，所以更适合第一步用“利率上升—行业收益反应”这个简单线性逻辑来解释。

但这里要注意：**如果论文后面要写得更严谨，不能只靠一个我挑出来的日期。** 单个事件容易被质疑“是不是事后挑了一个结果好看的日期”。所以更好的处理方式是：

1. 第一步用 2022-11-02 做教学式、第二讲对应的一元截面回归；
2. 稳健性检验换其他紧缩 FOMC 日期，比如 2022-06-15、2022-07-27、2022-09-21、2022-12-14；
3. 后面进阶版本把多个 FOMC 事件堆叠成“行业 × 事件”的面板或混合截面数据；
4. 最严格的版本不用“是否加息”定义冲击，而是用事件日前后 DGS2 或联邦基金期货的高频变化来度量货币政策 surprise。

所以 2022-11-02 的定位应该写成：**第一阶段的基准示范事件，而不是最终论文的唯一证据。**

[1]: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html "Kenneth R. French - Data Library"
[2]: https://fred.stlouisfed.org/series/DGS2 "Market Yield on U.S. Treasury Securities at 2-Year Constant Maturity, Quoted on an Investment Basis (DGS2) | FRED | St. Louis Fed"
[3]: https://www.federalreserve.gov/monetarypolicy/fomccalendars.htm "The Fed - Meeting calendars and information"
[4]: https://www.federalreserve.gov/newsevents/pressreleases/monetary20221102a.htm "Federal Reserve Board - Federal Reserve issues FOMC statement"


# 我问：这些只用链接就好了吗 C. FOMC事件日期与声明
2. 这些数据只用链接就好了吗？

按老师要求，压缩包里要包括原始数据；如果数据集较大，可以建一个 .txt 文件给出数据链接地址。这个要求在期末大作业说明里写得很明确：压缩包中须包括原始数据，若数据集较大可建立记事本给出数据链接地址。

我的建议是不要只放链接，最好这样放：

data/
  raw/
    49_Industry_Portfolios_daily_CSV.zip
    DGS2.csv
    fomc_events.csv
    data_sources.txt

  processed/
    step1_cross_section_49industries.csv

其中：

49_Industry_Portfolios_daily_CSV.zip：Kenneth French 原始行业收益数据；
DGS2.csv：FRED 下载下来的原始利率数据；
fomc_events.csv：我们自己整理的 FOMC 事件日期表；
data_sources.txt：写清楚所有数据网址；
processed/：放我们清洗后用于回归的数据。

这样老师看起来最舒服：既有原始数据，也有来源链接，也有可复现的处理结果。

3. FOMC 事件日期与声明应该怎么处理？

FOMC 事件不是一个像 DGS2 那样直接下载的数值序列。它更像是我们从美联储官网整理出来的“事件表”。

我们可以自己建一个 fomc_events.csv，第一步先这样：

event_date,policy_type,target_range_after,statement_url,notes
2022-11-02,tightening,3.75-4.00,https://www.federalreserve.gov/newsevents/pressreleases/monetary20221102a.htm,75bp rate hike

后面做多个事件时扩展成：

event_date,policy_type,target_range_after,statement_url,notes
2022-03-16,tightening,0.25-0.50,https://www.federalreserve.gov/newsevents/pressreleases/monetary20220316a.htm,25bp rate hike
2022-05-04,tightening,0.75-1.00,https://www.federalreserve.gov/newsevents/pressreleases/monetary20220504a.htm,50bp rate hike
2022-06-15,tightening,1.50-1.75,https://www.federalreserve.gov/newsevents/pressreleases/monetary20220615a.htm,75bp rate hike
2022-07-27,tightening,2.25-2.50,https://www.federalreserve.gov/newsevents/pressreleases/monetary20220727a.htm,75bp rate hike
2022-09-21,tightening,3.00-3.25,https://www.federalreserve.gov/newsevents/pressreleases/monetary20220921a.htm,75bp rate hike
2022-11-02,tightening,3.75-4.00,https://www.federalreserve.gov/newsevents/pressreleases/monetary20221102a.htm,75bp rate hike
2022-12-14,tightening,4.25-4.50,https://www.federalreserve.gov/newsevents/pressreleases/monetary20221214a.htm,50bp rate hike

然后 data_sources.txt 里写：

FOMC official calendars and statements:
https://www.federalreserve.gov/monetarypolicy/fomccalendars.htm

2022-11-02 FOMC statement:
https://www.federalreserve.gov/newsevents/pressreleases/monetary20221102a.htm


方法四：不用 FRED，改用美联储 H.15 官方数据

因为 DGS2 在 FRED 上的源头就是美联储 H.15 Selected Interest Rates，所以也可以直接从美联储下载。FRED 页面也标明该序列来源于 Board of Governors，并属于 H.15 Selected Interest Rates。

打开：

https://www.federalreserve.gov/datadownload/Choose.aspx?rel=H15

然后找：

Treasury Constant Maturities

下载 CSV All Observations。美联储数据下载页面中确实提供 Treasury Constant Maturities 的 CSV 全部观测下载。

这个文件里 2 年期的序列名称通常是：

RIFLGFCY02_N.B

也就是 2-year constant maturity Treasury yield。美联储 H.15 下载说明中显示它对应 “Market yield on U.S. Treasury securities at 2-year constant maturity, quoted on investment basis”。

我建议你现在采用的顺序

先手动用 FRED 页面下载；如果页面下载也不稳定，就用 R 的 quantmod::getSymbols("DGS2", src = "FRED")；如果 FRED 整体都打不开，再用美联储 H.15 数据下载页面。

在论文数据来源里，我们仍然可以写：

DGS2: FRED, Federal Reserve Bank of St. Louis.
Original source: Board of Governors of the Federal Reserve System, H.15 Selected Interest Rates.

这样来源是正规的。