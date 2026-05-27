# Processing NHANES 2003-2006 semi-raw accelerometer data (ActiGraph 7164, single-axis, 60-second epoch counts) using the same GGIR pipeline as the 2011-2014 cycles for cross-cycle comparability.

## Background

The NHANES 2003 and 2005 cycles released semi-raw accelerometer data nearly two decades ago:
- [2003 Cycle (PAXRAW_C)](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)
- [2005 Cycle (PAXRAW_D)](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.htm)

NCI has already released an [algorithm](https://epi.grants.cancer.gov/nhanes-pam/create.html) for processing this data, and there is a [great R package](https://github.com/vandomed/nhanesaccel) available as well. But my goal is different: I already processed the 2011-2014 cycles with GGIR, and I want results that are comparable across cycles using the same pipeline.

## Important Limitations

The 2003-2006 data were collected using the ActiGraph 7164, a single-axis accelerometer with 60-second epoch count data. Participants were instructed to remove the device while showering, swimming, or sleeping ([protocol](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.htm)).

As a result:
1. Most sleep detection and step count algorithms (e.g., Verisense) do not work on this data.
2. Circadian rhythm estimates are unreliable because wear time is generally insufficient. L5 is most likely to coincide with the non-wear period (during sleep), producing RA values near 1 in ~90% of cases, which is not meaningful. Several published studies have reported RA, M10, and L5 from the 2003-2006 NHANES cycles. These results should be interpreted with caution, as the hip-worn, wake-only wear protocol means L5 almost always falls in the non-wear (sleep) period, artificially inflating RA toward 1 regardless of the participant's actual rest-activity pattern.
3. In roughly 10% of participants who did not follow the removal protocol, circadian rhythm parameters may be estimable.

The output metric is ActiGraph activity counts (counts/min), not ENMO (mg), so this data is not directly comparable to the wrist-worn ENMO output of NHANES 2011-2014 or UK Biobank without a count-to-ENMO conversion.

## Step 0: Download

Download [PAXRAW_C](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/PAXRAW_C.zip) and [PAXRAW_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.zip).

## Step 1: NCI Calibration

Run the NCI calibration using [create.pam_perminute.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perminute.sas). NCI appears to have manually reviewed all 2003-2004 data and performed QC (see lines 122-156 of the SAS code).

This SAS file only covers 2003-2004. You can modify it for 2005-2006 by pointing the input to `paxraw_d` and the output to `pam_perminute_d`, leaving everything else unchanged: the manual QC block references only 2003-2004 SEQN (21xxx-30xxx), so it is a no-op on the 2005-2006 SEQN (31xxx-41xxx). After this step you will have two SAS files: `pam_perminute_c.sas7bdat` and `pam_perminute_d.sas7bdat`, totaling about 8.43 GB.

NCI's downstream code (see [create.pam_perday.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perday.sas)) gates on `PAXCAL=1` (monitor in calibration), not PAXSTAT. I replicate this by keeping only `PAXCAL==1` records in R (Step 2), which drops uncalibrated (PAXCAL=2) and unknown (PAXCAL=9) monitors.

## Step 2: Merge with Demographics

Merge age information from the DEMO files, since NCI uses age-specific cutoffs (see [create.pam_perday.sas](https://epi.grants.cancer.gov/nhanes-pam/create.pam_perday.sas), lines 206-250).

```r
pam_perminute_C <- read_sas("pam_perminute_c.sas7bdat")
pam_perminute_D <- read_sas("pam_perminute_d.sas7bdat")
demo_C <- read_xpt("DEMO_C.XPT")
demo_D <- read_xpt("DEMO_D.XPT")
demo <- bind_rows(demo_C, demo_D) %>% dplyr::select(SEQN, RIDAGEYR)
pam_perminute <- bind_rows(pam_perminute_C, pam_perminute_D %>% dplyr::select(-PAXSTEP))
pam_perminute <- inner_join(pam_perminute, demo, by = "SEQN")
pam_perminute <- filter(pam_perminute, PAXCAL == 1)   # NCI calibration gate
```

## Step 3: Generate ActiGraph-Style CSV Files

Format the data into ActiGraph GT3X-style CSV files using VM = count, matching the format used in the 2011-2014 pipeline. This ensures comparability and avoids writing a separate reader.

**Dating fix.** An earlier version mapped each minute's `PAXDAY` (day of week) directly to a fixed calendar date. This preserves the weekday but corrupts the recording order: a participant who started mid-week gets their days re-sorted by weekday rather than by wear sequence, so `setorder(date)` no longer reflects day 1 -> day 7. The fix derives the recording day from `PAXN` (the device's sequential minute counter, 1-10080, 1440 per day), and anchors day 1 to a calendar date whose weekday matches the real `PAXDAY`, then increments. This gives both correct chronological order and correct weekday. The real weekday is saved to an external lookup table keyed on (SEQN, day_in_rec).

**Timestamp format.** The data rows use ISO-8601 (`%Y-%m-%dT%H:%M:%SZ`), matching the [GGIR ActiGraph test file](https://github.com/wadpac/GGIRread/blob/main/inst/testfiles/ActiGraph13_timestamps_headers.csv). The header `Start Date` uses the ActiGraph convention `%m-%d-%Y`. GGIR (tested on v3.3.6) anchors time from the header `Start Date`, so `extEpochData_timeformat` in Step 4 must be set to `%m-%d-%Y %H:%M:%S` to match it; a slash format there raises a "time format does not match" error.

The code below processes adults (age >= 18). Age-specific cutoffs require separate processing for other age groups.

<details>
<summary>CSV generation code</summary>

```r
# PAXDAY (1=Sunday...7=Saturday) -> a 2000 calendar date with matching weekday
# 2000-01-02 is a Sunday
paxday_to_date <- setNames(
  as.Date(c("2000-01-02","2000-01-03","2000-01-04","2000-01-05",
            "2000-01-06","2000-01-07","2000-01-08")),
  as.character(1:7))
stopifnot(all(data.table::wday(paxday_to_date) == 1:7))  # 1=Sun..7=Sat

pam_perminuteWorking <- filter(pam_perminute, RIDAGEYR >= 18) %>%
  dplyr::select(-RIDAGEYR)
allSEQN <- unique(pam_perminuteWorking$SEQN)

weekday_records <- list()

for (seq in allSEQN) {
  temp <- filter(pam_perminuteWorking, SEQN == seq) %>% arrange(PAXN)

  # recording day from PAXN; anchor day 1 to its real weekday, then increment
  day_in_rec <- ((temp$PAXN - 1L) %/% 1440L) + 1L
  day1_date  <- paxday_to_date[as.character(temp$PAXDAY[which.min(temp$PAXN)])]

  out <- temp %>%
    mutate(
      day_in_rec = day_in_rec,
      Date = day1_date + (day_in_rec - 1L),
      TimeStamp = as.POSIXct(
        paste0(Date, " ", sprintf("%02d:%02d:00", PAXHOUR, PAXMINUT)),
        tz = "UTC"
      ),
      TimeStamp = format(TimeStamp, "%Y-%m-%dT%H:%M:%SZ"),
      vm    = ifelse(is.na(PAXINTEN), 0L, as.integer(PAXINTEN)),  # NA -> 0, GGIR nonwear catches it
      axis1 = vm, axis2 = 0, axis3 = 0  # dummy values; GGIR uses VM only
    )

  # save real weekday for downstream merge
  weekday_records[[as.character(seq)]] <-
    distinct(out, day_in_rec, weekday = temp$PAXDAY)

  body <- out %>% select(TimeStamp, axis1, axis2, axis3, vm)

  first_ts <- as.POSIXct(out$TimeStamp[1], format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  hdr <- c(
    "------------ Data Table File Created By Actigraph Link ActiLife v6.11.9 date format dd/MM/yyyy Filter Normal -----------,,,,,",
    "Serial Number: TAS1D48140206,,,,,",
    paste0("Start Time ", format(first_ts, "%H:%M:%S"), ",,,,,"),
    paste0("Start Date ", format(first_ts, "%m-%d-%Y"), ",,,,,"),
    "Epoch Period (hh:mm:ss) 00:00:60,,,,,",
    "Download Time 15:29:44,,,,,",
    "Download Date 09-19-2017,,,,,",
    "Current Memory Address: 0,,,,,",
    "Current Battery Voltage: 3.85     Mode = 13,,,,,",
    "--------------------------------------------------,,,,,"
  )

  file_path <- file.path("/NHANESGGIR/NHANES03040506Adult/", paste0(seq, ".csv"))
  writeLines(hdr, file_path)
  write.table(body, file_path, sep = ",", row.names = FALSE,
              col.names = TRUE, quote = FALSE, append = TRUE)
}

weekday_df <- bind_rows(weekday_records, .id = "SEQN")
arrow::write_parquet(weekday_df, "weekday_lookup.parquet")
```

</details>

## Step 4: Process with GGIR

Key design decisions for using GGIR on 2003-2006 data:

- **No calibration or imputation**: not feasible with single-axis count data (`do.cal = FALSE`, `do.imp = FALSE`)
- **Window sizes**: `c(60, 3600, 3600)` to match the NCI algorithm (60s native epoch, 1-hour non-wear and assessment windows)
- **Non-wear as sleep proxy**: uses GGIR's NotWorn algorithm for both `HASPT.algo` and `HASIB.algo`, since actual sleep detection is not possible with this data
- **10-hour minimum wear time** per valid day, matching NCI
- **Count thresholds are per-minute**, so no epoch rescaling is needed at 60s epochs (unlike the 5s GGIR vignette example)
- **NCI activity thresholds**: sedentary < 100 cpm, with NCI's age-specific moderate and vigorous cutoffs (Troiano 2008 for adults, Evenson 2008 for youth, taken directly from create.pam_perday.sas)
- **`extEpochData_timeformat = "%m-%d-%Y %H:%M:%S"`**: GGIR anchors time from the header `Start Date`, which uses this hyphenated format
- **NeishabouriCount columns are not usable**: they contain dummy values from the CSV formatting step; GGIR imports the VM column via `acc.metric`

Verify on a single SEQN that the first epoch in the `meta/basic` output lands on the expected calendar date (e.g. SEQN 21005, day 1 = Sunday, should start on 2000-01-02) and that ACC values are real counts before running the full cohort.

<details>
<summary>GGIR configuration (adults, age >= 18)</summary>

```r
GGIR(
  idloc                    = 6,
  datadir                  = "/NHANESGGIR/NHANES03040506Adult",
  outputdir                = "/NHANESGGIR/ResultNHANES03040506Adult",
  do.cal                   = FALSE,
  dataFormat               = "actigraph_csv",
  extEpochData_timeformat  = "%m-%d-%Y %H:%M:%S",
  windowsizes              = c(60, 3600, 3600),
  do.neishabouricounts     = FALSE,
  acc.metric               = "NeishabouriCount_vm",
  studyname                = "NHANES03040506",
  HASPT.algo               = "NotWorn",
  HASIB.algo               = "NotWorn",
  do.imp                   = FALSE,
  nonwear_range_threshold  = 100,
  HASPT.ignore.invalid     = NA,
  ignorenonwear            = FALSE,
  threshold.lig            = 100,
  threshold.mod            = 2020,   # NCI Troiano 2008, adult
  threshold.vig            = 5999,   # NCI Troiano 2008, adult
  includedaycrit           = 10,
  includenightcrit         = 10
)
```

</details>

<details>
<summary>Age-specific cutoffs (youth, age 6-17)</summary>

NCI uses age-specific moderate/vigorous count cutoffs for youth (Evenson 2008). Process each youth age in a separate GGIR run with its own `threshold.mod` / `threshold.vig`, keeping the rest of the configuration identical to the adult run.

```r
ages       <- 6:17
modThresh  <- c(1400,1515,1638,1770,1910,2059,2220,2393,2580,2781,3000,3239)  # NCI Evenson 2008
vigThresh  <- c(3758,3947,4147,4360,4588,4832,5094,5375,5679,6007,6363,6751)  # NCI Evenson 2008
names(modThresh) <- names(vigThresh) <- as.character(ages)

for (age in ages) {
  GGIR(
    idloc                    = 6,
    datadir                  = sprintf("/NHANESGGIR/NHANES03040506Age%d", age),
    outputdir                = sprintf("/NHANESGGIR/ResultNHANES03040506Age%d", age),
    do.cal                   = FALSE,
    dataFormat               = "actigraph_csv",
    extEpochData_timeformat  = "%m-%d-%Y %H:%M:%S",
    windowsizes              = c(60, 3600, 3600),
    do.neishabouricounts     = FALSE,
    acc.metric               = "NeishabouriCount_vm",
    studyname                = sprintf("NHANES03040506Age%d", age),
    HASPT.algo               = "NotWorn",
    HASIB.algo               = "NotWorn",
    do.imp                   = FALSE,
    nonwear_range_threshold  = 100,
    HASPT.ignore.invalid     = NA,
    ignorenonwear            = FALSE,
    threshold.lig            = 100,
    threshold.mod            = modThresh[as.character(age)],
    threshold.vig            = vigThresh[as.character(age)],
    includedaycrit           = 10,
    includenightcrit         = 10
  )
}
```

</details>

## Validation

I compared MVPA and sedentary time against the NCI SAS algorithm for participants aged >= 18 (raw correlation, no survey weights):

| Metric | Correlation with NCI |
|---|---|
| MVPA | 97.12% |
| Sedentary time | 85.43% |

MVPA correlation is high because at the 60-second epoch level, bout detection is unnecessary and the count threshold is straightforward. The lower sedentary correlation likely reflects differences in non-wear handling between GGIR and the NCI algorithm.

## Contact

If you have questions or would like access to the processed data (too large for GitHub), feel free to reach out on [LinkedIn](https://www.linkedin.com/in/xiaoyu-zong-0a733ba0/)
