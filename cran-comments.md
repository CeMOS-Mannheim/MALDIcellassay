## R CMD check results

0 errors | 0 warnings | 1 note

>Please reduce the length of the title to less than 65 characters.
>
>If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
>authors (year) <doi:...>
>authors (year, ISBN:...)
>or if those are not available: <https:...>
>with no space after 'doi:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")
>
>The Description field is intended to be a (one paragraph) description of what the package does and why it may be useful. Please add more details about the package functionality and implemented methods in your Description text.
>
>Please only write package and software names in single quotes.
>-> Omit them around: un-targeted
>
>Please fix and resubmit.
>
>Best,
>Benjamin Altmann

Dear CRAN-team
I took care about the comments above.

- The new title is: "Automated MALDI Cell Assays Using Dose-Response Curve Fitting" (63 character).
- DOI's have been added to the description field.
- The description field was extended: "Conduct automated cell-based assays using Matrix-Assisted Laser Desorption/Ionization (MALDI) methods for high-throughput screening of signals responsive to treatments. The package efficiently identifies high variance signals and fits dose-response curves to them. Quality metrics such as Z', V', log2FC, and CRS are provided for evaluating the potential of signals as biomarkers. The methodologies were introduced by Weigt et al. (2018) <doi:10.1038/s41598-018-29677-z> and refined by Unger et al. (2021) <doi:10.1038/s41596-021-00624-z>."
- Single quotes were removed from the title.

Best regards,
Thomas Enzlein

