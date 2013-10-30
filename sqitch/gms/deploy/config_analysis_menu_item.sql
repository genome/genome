-- Deploy config_analysis_menu_item
-- requires: config_schema
-- requires: config_set

BEGIN;

CREATE TABLE IF NOT EXISTS config.analysis_menu_item (
    id character varying(64) NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    name character varying(255) NOT NULL,
    configuration_set_id character varying(64) NOT NULL,
    CONSTRAINT analysis_menu_item_pk PRIMARY KEY (id),
    CONSTRAINT analysis_menu_config_set_fk FOREIGN KEY (configuration_set_id) REFERENCES config.set(id)
);


COMMIT;
