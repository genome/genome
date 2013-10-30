-- Deploy disk_file_summary
-- requires: disk_schema
-- requires: disk_allocation

BEGIN;

CREATE TABLE IF NOT EXISTS disk.file_summary (
    id character varying NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    created_by character varying NOT NULL,
    allocation_id character varying NOT NULL,
    file character varying NOT NULL,
    digest character varying,
    size_in_bytes bigint,
    is_symlink boolean NOT NULL,
    destination character varying,
    CONSTRAINT file_summary_pkey PRIMARY KEY (id),
    CONSTRAINT file_summary_allocation_id_file_key UNIQUE (allocation_id, file),
    CONSTRAINT file_summary_allocation_id_fkey FOREIGN KEY (allocation_id) REFERENCES disk.allocation(id)
);

COMMIT;
