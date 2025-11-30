# group_29_project

## Data source

The project used data from the article:

-   Frontiers in Immunology, “B Cell Biology,” vol. 10, Oct. 23, 2019. doi: 10.3389/fimmu.2019.02434

## Download data

The data was downloaded from the corresponding shiny app:

<http://bioinformin.cesnet.cz/CDmaps/>

**The workflow for the data retrieval was as follows:**

-   Go to the shiny app

-   Click on the `data` tab

-   Check the box *"Show individual samples data (by default only median values are displayed)"*

-   Then download the three files:

    -   `2. Blood: B- and T cells`

    -   `5. Thymus`

    -   `4. Tonsil: B-cell maturation`

-   As the names of downloaded files were automatically named after the current date, they were renamed to:

    -   `blood.csv`

    -   `thymus.csv`

    -   `tonsil.csv`

-   The three files were then uploaded to the folder `_raw` inside of the folder `data`
