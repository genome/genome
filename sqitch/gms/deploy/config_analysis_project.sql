-- Deploy config_analysis_project
-- requires: config_schema
-- requires: config_analysis_menu_item
-- requires: config_set
-- requires: model_model_group
-- requires: disk_allocation

BEGIN;

CREATE TABLE IF NOT EXISTS config.analysis_project (
    id character varying(64) NOT NULL,
    configuration_set_id character varying(64),
    created_by character varying(255) NOT NULL,
    analysis_menu_item_id character varying(64),
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    name character varying(255) NOT NULL,
    status character varying(255) NOT NULL,
    model_group_id character varying(64) NOT NULL,
    CONSTRAINT config_analysis_project_pk PRIMARY KEY(id),
    CONSTRAINT analysis_project_model_group FOREIGN KEY (model_group_id) REFERENCES model.model_group(id),
    CONSTRAINT config_analysis_project_analysis_menu_item_fk FOREIGN KEY (analysis_menu_item_id) REFERENCES config.analysis_menu_item(id),
    CONSTRAINT config_analysis_project_configruation_set_fk FOREIGN KEY (configuration_set_id) REFERENCES config.set(id)
);

COMMIT;
