-- Deploy config_set
-- requires: config_schema
-- requires: disk_allocation

BEGIN;

CREATE TABLE IF NOT EXISTS config.set (
    id character varying(64) NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    allocation_id character varying(64) NOT NULL,
    CONSTRAINT genome_config_set_pk PRIMARY KEY (id),
    CONSTRAINT config_set_allocation_fk FOREIGN KEY (allocation_id) REFERENCES disk.allocation(id)
);

COMMIT;
