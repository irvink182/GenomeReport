/* ==========================================================================
   GenomeReport
   Interactive report components
   ========================================================================== */


/* --------------------------------------------------------------------------
   Interactive tables
   -------------------------------------------------------------------------- */

document.querySelectorAll(".interactive-table").forEach((table) => {

    const container = table.closest(".table-container");

    if (!container) {
        return;
    }

    const tbody = table.querySelector("tbody");

    if (!tbody) {
        return;
    }

    const headers = Array.from(
        table.querySelectorAll("thead th")
    );

    const rows = Array.from(
        tbody.querySelectorAll("tr")
    );

    let filteredRows = [...rows];

    let currentPage = 1;
    let rowsPerPage = 25;

    let sortColumn = null;
    let sortAscending = true;


    /* ----------------------------------------------------------------------
       Controls
       ---------------------------------------------------------------------- */

    const controls = document.createElement("div");

    controls.className = "table-controls";


    /* Search */

    const searchGroup = document.createElement("div");

    searchGroup.className = "table-search";

    const searchLabel = document.createElement("label");

    searchLabel.textContent = "Search";

    const searchInput = document.createElement("input");

    searchInput.type = "search";
    searchInput.placeholder = "Search table...";

    searchLabel.appendChild(searchInput);
    searchGroup.appendChild(searchLabel);


    /* Rows per page */

    const pageSizeGroup = document.createElement("div");

    pageSizeGroup.className = "table-page-size";

    const pageSizeLabel = document.createElement("label");

    pageSizeLabel.textContent = "Rows";

    const pageSizeSelect = document.createElement("select");

    [10, 25, 50, 100].forEach((value) => {

        const option = document.createElement("option");

        option.value = value;
        option.textContent = value;

        if (value === rowsPerPage) {
            option.selected = true;
        }

        pageSizeSelect.appendChild(option);
    });

    pageSizeLabel.appendChild(pageSizeSelect);
    pageSizeGroup.appendChild(pageSizeLabel);


    controls.appendChild(searchGroup);
    controls.appendChild(pageSizeGroup);

    container.insertBefore(
        controls,
        table
    );


    /* ----------------------------------------------------------------------
       Footer
       ---------------------------------------------------------------------- */

    const footer = document.createElement("div");

    footer.className = "table-footer";


    const info = document.createElement("div");

    info.className = "table-info";


    const pagination = document.createElement("div");

    pagination.className = "table-pagination";


    footer.appendChild(info);
    footer.appendChild(pagination);

    container.appendChild(footer);


    /* ----------------------------------------------------------------------
       Value parsing
       ---------------------------------------------------------------------- */

    function parseValue(value) {

        const text = value.trim();

        if (
            text === ""
            || text.toLowerCase() === "na"
            || text.toLowerCase() === "nan"
        ) {
            return {
                type: "missing",
                value: null,
            };
        }

        const normalized = text
            .replace(/,/g, "")
            .replace(/%$/, "");

        const number = Number(normalized);

        if (!Number.isNaN(number)) {

            return {
                type: "number",
                value: number,
            };
        }

        return {
            type: "text",
            value: text.toLowerCase(),
        };
    }


    /* ----------------------------------------------------------------------
       Sorting
       ---------------------------------------------------------------------- */

    function sortRows() {

        if (sortColumn === null) {
            return;
        }

        filteredRows.sort((rowA, rowB) => {

            const valueA = parseValue(
                rowA.children[sortColumn].textContent
            );

            const valueB = parseValue(
                rowB.children[sortColumn].textContent
            );


            /* Missing values always last */

            if (
                valueA.type === "missing"
                && valueB.type === "missing"
            ) {
                return 0;
            }

            if (valueA.type === "missing") {
                return 1;
            }

            if (valueB.type === "missing") {
                return -1;
            }


            let comparison;


            /* Numerical comparison */

            if (
                valueA.type === "number"
                && valueB.type === "number"
            ) {

                comparison =
                    valueA.value - valueB.value;

            }

            /* Text comparison */

            else {

                comparison = String(
                    valueA.value
                ).localeCompare(
                    String(valueB.value),
                    undefined,
                    {
                        numeric: true,
                        sensitivity: "base",
                    }
                );
            }


            return (
                sortAscending
                ? comparison
                : -comparison
            );
        });
    }


    headers.forEach((header, index) => {

        header.classList.add(
            "sortable"
        );

        header.setAttribute(
            "tabindex",
            "0"
        );

        header.setAttribute(
            "role",
            "button"
        );


        function activateSort() {

            if (sortColumn === index) {

                sortAscending =
                    !sortAscending;

            } else {

                sortColumn = index;
                sortAscending = true;
            }


            headers.forEach((item) => {

                item.classList.remove(
                    "sort-asc",
                    "sort-desc"
                );

                item.removeAttribute(
                    "aria-sort"
                );
            });


            header.classList.add(
                sortAscending
                ? "sort-asc"
                : "sort-desc"
            );

            header.setAttribute(
                "aria-sort",
                sortAscending
                    ? "ascending"
                    : "descending"
            );


            currentPage = 1;

            renderTable();
        }


        header.addEventListener(
            "click",
            activateSort
        );


        header.addEventListener(
            "keydown",
            (event) => {

                if (
                    event.key === "Enter"
                    || event.key === " "
                ) {

                    event.preventDefault();

                    activateSort();
                }
            }
        );
    });


    /* ----------------------------------------------------------------------
       Search
       ---------------------------------------------------------------------- */

    function filterRows() {

        const query = searchInput.value
            .trim()
            .toLowerCase();

        if (query === "") {

            filteredRows = [...rows];

        } else {

            filteredRows = rows.filter((row) => {

                return row.textContent
                    .toLowerCase()
                    .includes(query);
            });
        }

        currentPage = 1;

        renderTable();
    }


    searchInput.addEventListener(
        "input",
        filterRows
    );


    /* ----------------------------------------------------------------------
       Rows per page
       ---------------------------------------------------------------------- */

    pageSizeSelect.addEventListener(
        "change",
        () => {

            rowsPerPage = Number(
                pageSizeSelect.value
            );

            currentPage = 1;

            renderTable();
        }
    );


    /* ----------------------------------------------------------------------
       Pagination
       ---------------------------------------------------------------------- */

    function createPageButton(
        label,
        page,
        disabled=false,
        active=false
    ) {

        const button = document.createElement(
            "button"
        );

        button.type = "button";
        button.textContent = label;

        button.disabled = disabled;

        if (active) {
            button.classList.add("active");
        }

        button.addEventListener(
            "click",
            () => {

                currentPage = page;

                renderTable();
            }
        );

        return button;
    }


    function renderPagination(
        totalPages
    ) {

        pagination.innerHTML = "";

        if (totalPages <= 1) {
            return;
        }


        /* Previous */

        pagination.appendChild(
            createPageButton(
                "Previous",
                currentPage - 1,
                currentPage === 1
            )
        );


        /*
         * Keep pagination compact.
         *
         * Example:
         *
         * 1  ...  6  7  8  ...  42
         */

        const pages = new Set([
            1,
            totalPages,
            currentPage - 1,
            currentPage,
            currentPage + 1,
        ]);

        const validPages = Array.from(pages)
            .filter(
                (page) =>
                    page >= 1
                    && page <= totalPages
            )
            .sort(
                (a, b) => a - b
            );


        let previousPage = null;

        validPages.forEach((page) => {

            if (
                previousPage !== null
                && page - previousPage > 1
            ) {

                const ellipsis =
                    document.createElement("span");

                ellipsis.className =
                    "pagination-ellipsis";

                ellipsis.textContent = "…";

                pagination.appendChild(
                    ellipsis
                );
            }


            pagination.appendChild(
                createPageButton(
                    String(page),
                    page,
                    false,
                    page === currentPage
                )
            );

            previousPage = page;
        });


        /* Next */

        pagination.appendChild(
            createPageButton(
                "Next",
                currentPage + 1,
                currentPage === totalPages
            )
        );
    }


    /* ----------------------------------------------------------------------
       Render table
       ---------------------------------------------------------------------- */

    function renderTable() {

        sortRows();


        const totalRows =
            filteredRows.length;

        const totalPages = Math.max(
            1,
            Math.ceil(
                totalRows / rowsPerPage
            )
        );


        if (currentPage > totalPages) {
            currentPage = totalPages;
        }


        const start =
            (currentPage - 1)
            * rowsPerPage;

        const end = Math.min(
            start + rowsPerPage,
            totalRows
        );


        /*
         * Detach all rows before adding
         * the rows for the current page.
         */

        tbody.innerHTML = "";


        filteredRows
            .slice(start, end)
            .forEach((row) => {

                tbody.appendChild(row);
            });


        /* Information */

        if (totalRows === 0) {

            info.textContent =
                "No matching entries";

        } else {

            info.textContent =
                `Showing ${start + 1}–${end} `
                + `of ${totalRows} entries`;

            if (totalRows !== rows.length) {

                info.textContent +=
                    ` (filtered from `
                    + `${rows.length})`;
            }
        }


        renderPagination(
            totalPages
        );
    }


    /* Initial rendering */

    renderTable();
});


/* --------------------------------------------------------------------------
   Figure enlargement
   -------------------------------------------------------------------------- */

const figureImages = document.querySelectorAll(
    ".figure-container img"
);

if (figureImages.length > 0) {

    const modal = document.createElement(
        "div"
    );

    modal.className = "figure-modal";

    modal.setAttribute(
        "aria-hidden",
        "true"
    );


    const modalContent =
        document.createElement("div");

    modalContent.className =
        "figure-modal-content";


    const closeButton =
        document.createElement("button");

    closeButton.type = "button";

    closeButton.className =
        "figure-modal-close";

    closeButton.setAttribute(
        "aria-label",
        "Close enlarged figure"
    );

    closeButton.textContent = "×";


    const modalImage =
        document.createElement("img");

    modalImage.alt = "";


    modalContent.appendChild(
        closeButton
    );

    modalContent.appendChild(
        modalImage
    );

    modal.appendChild(
        modalContent
    );

    document.body.appendChild(
        modal
    );


    function closeModal() {

        modal.classList.remove(
            "open"
        );

        modal.setAttribute(
            "aria-hidden",
            "true"
        );

        modalImage.src = "";
        modalImage.alt = "";
    }


    figureImages.forEach((image) => {

        image.classList.add(
            "expandable-figure"
        );

        image.setAttribute(
            "tabindex",
            "0"
        );

        image.setAttribute(
            "role",
            "button"
        );


        function openFigure() {

            modalImage.src =
                image.src;

            modalImage.alt =
                image.alt || "Enlarged figure";

            modal.classList.add(
                "open"
            );

            modal.setAttribute(
                "aria-hidden",
                "false"
            );

            closeButton.focus();
        }


        image.addEventListener(
            "click",
            openFigure
        );


        image.addEventListener(
            "keydown",
            (event) => {

                if (
                    event.key === "Enter"
                    || event.key === " "
                ) {

                    event.preventDefault();

                    openFigure();
                }
            }
        );
    });


    closeButton.addEventListener(
        "click",
        closeModal
    );


    modal.addEventListener(
        "click",
        (event) => {

            if (event.target === modal) {
                closeModal();
            }
        }
    );


    document.addEventListener(
        "keydown",
        (event) => {

            if (
                event.key === "Escape"
                && modal.classList.contains(
                    "open"
                )
            ) {

                closeModal();
            }
        }
    );
}